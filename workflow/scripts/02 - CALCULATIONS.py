__author__ = "Graeme Ford"
__credits__ = ["Graeme Ford", "Prof. Michael S. Pepper", "Prof. Fourie Joubert", "Antionette Colic"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "rob@spot.colorado.edu"
__status__ = "Beta"



# %%
# Import dependancies
import gzip
import io
import json
import sys
import time
from os.path import join
from typing import Generator, Tuple

import pandas as pd
import requests
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

pd.options.mode.chained_assignment = None

# %%

# Set constants and functions to be used:
# locations = snakemake.config['locations'].keys()
with open(join("..", "..", "config", "config.json")) as f:
    config = json.load(f)

clusters = config["cluster"]["clusters"]
genes = config["locations"]
geneSummary = dict()
for cluster in config["cluster"]["clusters"]:
    geneSummary[cluster] = pd.read_excel(
        join("..", "..", "config", "{}".format(config["cluster"]["file"]))
    )[["ID", cluster]]
popsFile = pd.read_excel(join("..", "..", "config", "Clusters.xlsx"))
populations = ["AFR", "AMR", "EUR", "EAS", "SAS"]

refPop = "AFR"
compPop = ["AMR", "EUR", "EAS", "SAS"]


#  %%

#  Set POST Variables and Headers:
locations = config["locations"].keys()
populations = ["AFR", "AMR", "EUR", "EAS", "SAS"]
endpoint = "https://rest.ensembl.org/vep/homo_sapiens/region/"
headers = {"Content-Type": "application/json", "Accept": "application/json"}


def generate_params(key):
    params = {
        "hgvs": True,
        "CADD": True,
        "Phenotypes": True,
        "domains": True,
        # 'canonical': True,
        "refseq": True,
        "LoF": True,
        "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred",
        # 'transcript_id': config['locations'][key]["GRCh38"]["transcript_id"]
    }
    return params
    # if key == 'CYP2A6':
    #     return params | dict(transcrcipt_id="NM_000762.6")
    # if key == 'CYP2B6':
    #     return params | dict(transcript_id="NM_000767.5")
    # if key == 'UGT2B7':
    #     return params | dict(transcript_id="NM_001074.4")


def freq(alt: int, total: int) -> int:
    """Calculate percentage frequency for alleles based on count data.
    Args:
        alt (int): Number of observations of the alternate allele.
        total (int): Total number of observations of any allele.
    Returns:
        int: decimal percentage allele frequency.
    """
    if alt != 0 and total != 0:
        return row["ALT_CTS"] / row["OBS_CT"]
    if (alt == 0 and total != 0) or (alt == 0 and total == 0):
        return 0
    if alt != 0 and total == 0:
        return 999


condelTests = ["SIFT", "PolyPhen"]
cutoff = {"SIFT": 0.15, "PolyPhen": 0.28, "Condel": 0.46}
maximums = {"SIFT": 1, "PolyPhen": 1}
probabilities = {
    "SIFT": pd.read_csv(
        join("..", "..", "config", "CONDEL", "sift.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
    "PolyPhen": pd.read_csv(
        join("..", "..", "config", "CONDEL", "polyphen.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
}


# %%

# Define formulas for later use:
def read_vcf(path: str) -> pd.DataFrame:
    """Function to import a VCF file as a Pandas dataframe. Does not include Genotype columns. CREDIT: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    Args:
        path (str): Path to the VCF file in question
    Returns:
        Dataframe: A Pandas Dataframe of the VCF files content.
    """
    with gzip.open(path, "rt") as f:
        lines = [l for l in f if not l.startswith("##")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})


def generate_notation(row: pd.Series, gene: str) -> str:
    """Parses a row and returns the HGVS notation to query E! Ensembl.
    Args:
        row (Series): A row (generated from .iterrows() method) for which notation is needed.
    Returns:
        str: HGVS notation for the given variant row.
    """
    alleles = row["ALT"].split(",")
    notation = list(str())
    for allele in alleles:
        if (len(row["REF"]) > len(allele)) | (len(row["REF"]) == len(allele)):
            stop_coordinates = int(row["POS"]) + (len(row["REF"]) - 1)
            notation.append(
                "{S}:{P}-{P2}:1/{I}".format(
                    S=row["CHROM"], P=row["POS"], P2=stop_coordinates, I=allele
                )
            )
        elif len(row["REF"]) < len(allele):
            stop_coordinates = int(row["POS"]) + len(row["REF"])
            notation.append(
                "{S}:{P}-{P2}:{ST}/{I}".format(
                    S=row["CHROM"],
                    P=row["POS"],
                    P2=stop_coordinates,
                    I=allele,
                    ST=config["locations"][gene]["GRCh38"]["strand"],
                )
            )
        return notation


def chunk(dataset: pd.DataFrame, size: int) -> Generator:
    """Renders a generator to yield chunks of data from the original file.
    Args:
        dataset (pd.DataFrame): The dataset to divide into chunks.
        size (int): The maximum size of the resulting chunks.
    Returns:
        [type]: [description]
    Yields:
        Generator: [description]
    """
    return (dataset[pos : pos + size] for pos in range(0, len(dataset), size))


def merge(value: set) -> str:
    """Generate row value to save in df
    Args:
        value (set): The set you wish to convert
    Returns:
        str: The converted value, ready to store.
    """
    r = list(value)
    if len(r) == 1:
        return str(r[0])
    else:
        return " | ".join(value)


def CONDEL(sift: int, polyphen: int) -> Tuple[int, str]:
    """Function to calculate a weighted average score accross SIFT and PolyPhenV2 scores.
    Args:
        sift (int): A number indicating the likelyhood a variant affects protein function (Based on sequence homology and protein function)
        polyphen (int): A number indicating impact of a variant on protein structure using straightforward physical comparative methods
    Returns:
        Tuple[int, str]: A weighted score which favours variants with scores further away from each set cutoff (i.e. less ambiguity regarding the prediction) as well as a string, indicating the cutoff verdict
    """

    score = int()

    if sift <= cutoff["SIFT"]:
        score += (1 - sift / maximums["SIFT"]) * (
            1
            - probabilities["SIFT"]
            .loc[probabilities["SIFT"]["Score"] == sift, "Normal"]
            .values[0]
        )
    else:
        score += (1 - sift / maximums["SIFT"]) * (
            1
            - probabilities["SIFT"]
            .loc[probabilities["SIFT"]["Score"] == sift, "Deleterious"]
            .values[0]
        )

    if polyphen >= cutoff["PolyPhen"]:
        score += (polyphen / maximums["PolyPhen"]) * (
            1
            - probabilities["PolyPhen"]
            .loc[probabilities["PolyPhen"]["Score"] == polyphen, "Normal"]
            .values[0]
        )
    else:
        score += (polyphen / maximums["PolyPhen"]) * (
            1
            - probabilities["PolyPhen"]
            .loc[probabilities["PolyPhen"]["Score"] == polyphen, "Deleterious"]
            .values[0]
        )

    # This acocunts for number of VEP teests used. Since we are only using 2, we use 2...
    score = score / 2

    if score >= 0.469:
        pred = "deleterious"
    elif score >= 0 and score < 0.469:
        pred = "neutral"

    return (score, pred)


def Fishers(input: dict, refPop: str, compPop: list):
    """Runs a row-wise Fisher's Exact Test between the two listed populations
    Parameters:
    input (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data. 
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """
    for key, dataset in input.items():
        dataset["OR"] = dataset["Count"][["ID", "POS", "REF", "ALT"]]
        dataset["P"] = dataset["Count"][["ID", "POS", "REF", "ALT"]]
        for pop in compPop:
            for key2, direction in {
                "L": "less",
                "G": "greater",
                "T": "two-sided",
            }.items():
                oLabel = "{refPop}_{tail}_{pop}".format(
                    refPop=refPop, tail=key2, pop=pop
                )
                pLabel = "{refPop}_{tail}_{pop}".format(
                    refPop=refPop, tail=key2, pop=pop
                )
                dataset["OR"][oLabel] = None
                dataset["P"][pLabel] = None
                for index, row in dataset["Count"].iterrows():
                    contingency = [
                        [
                            row["{}_ac".format(refPop)],
                            row["{}_tc".format(refPop)] - row["{}_ac".format(refPop)],
                        ],
                        [
                            row["{}_ac".format(pop)],
                            row["{}_tc".format(pop)] - row["{}_ac".format(pop)],
                        ],
                    ]
                    oVal, pVal = fisher_exact(contingency, alternative=direction)
                    dataset["P"].loc[index, pLabel] = pVal
                    dataset["OR"].loc[index, oLabel] = oVal
                dataset["P"][pLabel] = multipletests(
                    dataset["P"][pLabel], method="bonferroni"
                )[1]
        columnsToDrop = list()
        # for pop in populations:
        #     columnsToDrop.append("{}_ac".format(pop))
        #     columnsToDrop.append("{}_tc".format(pop))
        # dataset.drop(columns=columnsToDrop ,inplace=True)


# %%

#  Import Data:
data = dict()
for location in locations:
    data[location] = read_vcf(
        join("..", "..", "results", "ALL_{}.vcf.gz".format(location))
    )

# %%

# Sub-Divide data:
data_generator = dict()
for dataset in data:
    data_generator[dataset] = chunk(data[dataset], 100)

# %%

# Compile and format request bodies:
data_to_send = dict()

# Iterate through each gene:
for dataset in data_generator:
    data_to_send[dataset] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in data_generator[dataset]:
        temp_list = list()

        # Iterate through each row in the chunk and add the HGVS notation to the list:
        for index, row in chunk.iterrows():
            temp_list.extend(generate_notation(row, dataset))
        data_to_send[dataset].append(dict(variants=temp_list))

# %%

# Perform API calls:
data_received = dict()
for dataset_key, dataset in data_to_send.items():
    data_received[dataset_key] = list()
    for index, chunk in enumerate(dataset):
        requesting = True
        temp_list = list()
        while requesting:
            r = requests.post(
                endpoint,
                headers=headers,
                data=json.dumps(chunk),
                params=generate_params(dataset_key),
            )
            if not r.ok:
                time.sleep(2)
            else:
                requesting = False
                decoded = r.json()
                data_received[dataset_key] = data_received[dataset_key] + decoded

# %%

# Iterate through each response and compile its values:
supplementary = dict()
for dataset_key, dataset in data_received.items():
    supplementary[dataset_key] = data[dataset_key][["ID", "POS", "REF", "ALT"]]
    new_columns = [
        "Co-Located Variant",
        "Transcript ID",
        "Transcript Strand",
        "Existing Variation",
        "Start Coordinates",
        "Consequence",
        "Diseases",
        "Biotype",
        "CADD_PHRED",
        # 'LoFtool',
        "input",
        "SIFT_score",
        "SIFT_pred",
        "Polyphen_score",
        "Polyphen_pred",
        "CONDEL",
        "CONDEL_pred",
    ]
    for column in new_columns:
        supplementary[dataset_key][column] = None

    for key, chunk in data_received.items():
        for variant in chunk:
            # row = supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] == int(variant['start'])]
            supplementary[dataset_key].loc[
                supplementary[dataset_key]["POS"] == int(variant["start"]),
                "Start Coordinates",
            ] = variant["start"]
            supplementary[dataset_key].loc[
                supplementary[dataset_key]["POS"] == int(variant["start"]), "input"
            ] = variant["input"]

            co_variants = list()
            if "colocated_variants" in variant:
                supplementary[dataset_key].loc[
                    supplementary[dataset_key]["POS"] == int(variant["start"]),
                    "Co-Located Variant",
                ] = True
                for colocated_variant in variant["colocated_variants"]:
                    co_variants.append(colocated_variant["id"])
            else:
                supplementary[dataset_key].loc[
                    supplementary[dataset_key]["POS"] == int(variant["start"]),
                    "Co-Located Variant",
                ] = False
                co_variants.append("-")
            supplementary[dataset_key].loc[
                supplementary[dataset_key]["POS"] == int(variant["start"]),
                "Existing Variation",
            ] = "| ".join(co_variants)

            if "transcript_consequences" in variant:
                for transcript in config["locations"][dataset_key]["GRCh38"][
                    "transcript_id"
                ]:
                    consequence = next(
                        (
                            n
                            for n in variant["transcript_consequences"]
                            if n["transcript_id"] == transcript
                        ),
                        None,
                    )
                    if consequence is not None:
                        phenotype = set()

                        if "phenotypes" in consequence:
                            for instance in consequence["phenotypes"]:
                                phenotype.add(instance["phenotype"])

                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "SIFT_score",
                        ] = (
                            consequence["sift_score"]
                            if ("sift_score" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "SIFT_pred",
                        ] = (
                            consequence["sift_prediction"]
                            if ("sift_prediction" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Polyphen_score",
                        ] = (
                            str(consequence["polyphen_score"])
                            if ("polyphen_score" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Polyphen_pred",
                        ] = (
                            consequence["polyphen_prediction"]
                            if ("polyphen_prediction" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Diseases",
                        ] = merge(phenotype)
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Consequence",
                        ] = merge(consequence["consequence_terms"])
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Transcript ID",
                        ] = consequence["transcript_id"]
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Biotype",
                        ] = (
                            consequence["biotype"]
                            if ("biotype" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "CADD_PHRED",
                        ] = (
                            consequence["cadd_phred"]
                            if ("cadd_phred" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Transcript Strand",
                        ] = (
                            consequence["strand"] if ("strand" in consequence) else None
                        )
                        # row['LoFtool'] = consequence['loftool'] if 'loftool' in consequence else None
                        if (
                            "sift_score" in consequence
                            and "polyphen_score" in consequence
                        ):
                            s, p = CONDEL(
                                int(consequence["sift_score"]),
                                consequence["polyphen_score"],
                            )
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL",
                            ] = s
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL_pred",
                            ] = p
                        else:
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL",
                            ] = None
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL_pred",
                            ] = None
                        break

                    else:
                        consequence = variant["transcript_consequences"][0]
                        phenotype = set()

                        if "phenotypes" in consequence:
                            for instance in consequence["phenotypes"]:
                                phenotype.add(instance["phenotype"])

                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "SIFT_score",
                        ] = (
                            consequence["sift_score"]
                            if ("sift_score" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "SIFT_pred",
                        ] = (
                            consequence["sift_prediction"]
                            if ("sift_prediction" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Polyphen_score",
                        ] = (
                            str(consequence["polyphen_score"])
                            if ("polyphen_score" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Polyphen_pred",
                        ] = (
                            consequence["polyphen_prediction"]
                            if ("polyphen_prediction" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Diseases",
                        ] = merge(phenotype)
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Consequence",
                        ] = merge(consequence["consequence_terms"])
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Transcript ID",
                        ] = consequence["transcript_id"]
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Biotype",
                        ] = (
                            consequence["biotype"]
                            if ("biotype" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "CADD_PHRED",
                        ] = (
                            consequence["cadd_phred"]
                            if ("cadd_phred" in consequence)
                            else None
                        )
                        supplementary[dataset_key].loc[
                            supplementary[dataset_key]["POS"] == int(variant["start"]),
                            "Transcript Strand",
                        ] = (
                            consequence["strand"] if ("strand" in consequence) else None
                        )
                        # row['LoFtool'] = consequence['loftool'] if 'loftool' in consequence else None
                        if (
                            "sift_score" in consequence
                            and "polyphen_score" in consequence
                        ):
                            s, p = CONDEL(
                                int(consequence["sift_score"]),
                                consequence["polyphen_score"],
                            )
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL",
                            ] = s
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL_pred",
                            ] = p
                        else:
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL",
                            ] = None
                            supplementary[dataset_key].loc[
                                supplementary[dataset_key]["POS"]
                                == int(variant["start"]),
                                "CONDEL_pred",
                            ] = None
                        break

# %%

# Save formatted results to CSV
for cluster in clusters:
    for gene in locations:
        supplementary[gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                "{}_VEP.csv".format(gene),
            ),
            sep="\t",
            index=False,
        )
        # supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)

# %%

# FREQUENCY CALCULATIONS
# %%

# Retrive VEP data:
supplementary = dict()

for cluster in clusters:
    supplementary[cluster] = dict()
    for gene in genes:
        supplementary[cluster][gene] = pd.read_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                "{}_VEP.csv".format(gene),
            ),
            sep="\t",
        )[["ID", "POS", "REF", "ALT"]]

# %%


# Calculate frequencies
frequency_data = dict()
fishers_data = dict()


for cluster in clusters:
    frequency_data[cluster] = dict()
    fishers_data[cluster] = dict()
    for gene in genes:
        fishers_data[cluster][gene] = supplementary[cluster][gene][
            ["ID", "POS", "REF", "ALT"]
        ]
        for pop in popsFile[cluster].unique():
            try:
                # supplementary[cluster][gene][pop] = 0
                fishers_data[cluster][gene]["{}_ac".format(pop)] = 0
                fishers_data[cluster][gene]["{}_tc".format(pop)] = 0
                frequency_data[cluster][gene] = pd.read_csv(
                    join(
                        "..",
                        "..",
                        "results",
                        cluster,
                        "ALL_{gene}.{pop}.acount".format(gene=gene, pop=pop),
                    ),
                    delimiter="\t",
                ).rename(columns={"#CHROM": "CHROM"})
                for index, row in frequency_data[cluster][gene].iterrows():
                    supplementary[cluster][gene].loc[
                        supplementary[cluster][gene]["ID"] == row["ID"], pop
                    ] = freq(row["ALT_CTS"], row["OBS_CT"])
                    fishers_data[cluster][gene].loc[
                        supplementary[cluster][gene]["ID"] == row["ID"],
                        "{pop}_ac".format(pop=pop),
                    ] = row["ALT_CTS"]
                    fishers_data[cluster][gene].loc[
                        supplementary[cluster][gene]["ID"] == row["ID"],
                        "{pop}_tc".format(pop=pop),
                    ] = row["OBS_CT"]
            except:
                pass

# %%

# Save the resulting dataframe back to its excel file:
for cluster in clusters:
    for gene in genes:
        supplementary[cluster][gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                "{}_Freq.csv".format(gene),
            ),
            index=False,
            sep="\t",
        )
        fishers_data[cluster][gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                "{}_Count.csv".format(gene),
            ),
            index=False,
            sep="\t",
        )


# %%

# FISHERS EXACT

# %%

# Load the Supplementary Table
supplementary = dict()

for cluster in clusters:
    supplementary[cluster] = dict()
    for gene in genes:
        supplementary[cluster][gene] = dict()
        supplementary[cluster][gene]["Count"] = pd.read_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                "{}_Count.csv".format(gene),
            ),
            sep="\t",
        )
        supplementary[cluster][gene]["OR"] = dict()
        supplementary[cluster][gene]["P"] = dict()

# %%

# Run Fisher's Exact test
Fishers(supplementary["SUPER"], refPop, compPop)

# %%

# Save Fishers Data to CSV
for gene in genes:
    supplementary["SUPER"][gene]["P"].to_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            "SUPER",
            "{}_FishersP.csv".format(gene),
        ),
        sep="\t",
        index=False,
    )
    supplementary["SUPER"][gene]["OR"].to_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            "SUPER",
            "{}_FishersOR.csv".format(gene),
        ),
        sep="\t",
        index=False,
    )
# %%

