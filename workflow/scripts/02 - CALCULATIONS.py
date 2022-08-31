#!/usr/bin/env python
"""
A Python script designed to run Frequency, Fishers Exact Test and
Variant Effect Prediction calculations and calls.
"""
# %%
# Import dependancies
import gzip
import io
import json
import time
from os import makedirs
from os.path import exists, join
from typing import Generator, Tuple

import pandas as pd
import requests
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

pd.options.mode.chained_assignment = None


__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

# %%

####################################################
############ IMPORT ALL CONTEXTUAL DATA ############
####################################################

# Import all clustering provided for this analysis from ../../input/samples.csv
CLUSTER_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))
CLUSTERS = set(CLUSTER_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")

LOCATIONS = pd.read_csv(join("..", "..", "input", "locations.csv"))
LOCATION_NAMES = LOCATIONS["location_name"].unique().tolist()

TRANSCRIPTS = pd.read_csv(join("..", "..", "input", "transcripts.csv"))

# %%
# geneSummary = dict()
# for cluster in clusters:
#     geneSummary[cluster] = pd.read_excel(
#         join("..", "..", "config", "{}".format(config["cluster"]["file"]))
#     )[["ID", cluster]]
# popsFile = pd.read_excel(join("..", "..", "config", "Clusters.xlsx"))
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]
ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# %%
def generate_params(transcription_paramaters: DataFrame, input_gene: str) -> dict:
    """Returns generated parameters suitable for the E! Ensembl VEP tool API.

    Returns:
        dict: A dictionary of parameters suitable for the VEP API.
    """
    params = {
        "hgvs": True,
        "CADD": True,
        "Phenotypes": True,
        "domains": True,
        # 'canonical': True,
        "refseq": True,
        "LoF": True,
        "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred",
        "transcript_id": TRANSCRIPTS.loc[
            TRANSCRIPTS["gene_name"] == input_gene, "transcript_id"
        ].tolist(),
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


# %%
CONDEL_TESTS = ["SIFT", "PolyPhen"]
TEST_MINIMUM_CUTOFFS = {"SIFT": 0.15, "PolyPhen": 0.28, "Condel": 0.46}
TEST_MAXIMUM_CUTOFFS = {"SIFT": 1, "PolyPhen": 1}
PROBABILITIES = {
    "SIFT": pd.read_csv(
        join("..", "..", "resources", "CONDEL", "sift.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
    "PolyPhen": pd.read_csv(
        join("..", "..", "resources", "CONDEL", "polyphen.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
}


def directory_exists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not exists(path):
        makedirs(path)


# Define formulas for later use:
def read_vcf(path: str) -> pd.DataFrame:
    """
    Function to import a VCF file as a Pandas dataframe. Does not include Genotype columns.
    CREDIT: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

    Args:
        path (str): Path to the VCF file in question
    Returns:
        Dataframe: A Pandas Dataframe of the VCF files content.
    """
    with gzip.open(path, "rt") as file:
        lines = [l for l in file if not l.startswith("##")]
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


def generate_notation(input_row: pd.Series, input_dataset: str) -> str:
    """
    Parses a row and returns the HGVS notation to query E! Ensembl.

    Args:
        row (Series): A row (generated from .iterrows() method) for which notation is needed.
    Returns:
        str: HGVS notation for the given variant row.
    """
    alleles = input_row["ALT"].split(",")
    notation = list()
    for allele in alleles:
        if (len(input_row["REF"]) > len(allele)) | (
            len(input_row["REF"]) == len(allele)
        ):
            stop_coordinates = int(input_row["POS"]) + (len(input_row["REF"]) - 1)
            notation.append(
                f"{input_row['CHROM']}:{input_row['POS']}-{stop_coordinates}:1/{allele}"
            )
        elif len(input_row["REF"]) < len(allele):
            stop_coordinates = int(input_row["POS"]) + len(input_row["REF"])
            notation.append(
                f"{input_row['CHROM']}:{input_row['POS']}-{stop_coordinates}:{LOCATIONS.loc[LOCATIONS['location_name'] == input_dataset, 'strand'].item()}/{allele}"
            )
        return notation


def chunk(input_data: pd.DataFrame, size: int) -> Generator:
    """Renders a generator to yield chunks of data from the original file.
    Args:
        dataset (pd.DataFrame): The dataset to divide into chunks.
        size (int): The maximum size of the resulting chunks.
    Returns:
        [type]: [description]
    Yields:
        Generator: [description]
    """
    return (input_data[pos : pos + size] for pos in range(0, len(input_data), size))


def merge(input_row: set) -> str:
    """Generate row value to save in df
    Args:
        value (set): The set you wish to convert
    Returns:
        str: The converted value, ready to store.
    """
    row_value_as_list = list(input_row)
    if len(row_value_as_list) == 1:
        return str(row_value_as_list[0])
    else:
        return " | ".join(input_row)


def condel_weighted_score(sift: int, polyphen: int) -> Tuple[int, str]:
    """Function to calculate a weighted average score accross SIFT and PolyPhenV2 scores.
    Args:
        sift (int): A number indicating the likelyhood a variant affects protein function (Based on sequence homology and protein function)
        polyphen (int): A number indicating impact of a variant on protein structure using straightforward physical comparative methods
    Returns:
        Tuple[int, str]: A weighted score which favours variants with scores further away from each set cutoff (i.e. less ambiguity regarding the prediction) as well as a string, indicating the cutoff verdict
    """

    score = int()

    if sift <= TEST_MINIMUM_CUTOFFS["SIFT"]:
        score += (1 - sift / TEST_MAXIMUM_CUTOFFS["SIFT"]) * (
            1
            - PROBABILITIES["SIFT"]
            .loc[PROBABILITIES["SIFT"]["Score"] == sift, "Normal"]
            .values[0]
        )
    else:
        score += (1 - sift / TEST_MAXIMUM_CUTOFFS["SIFT"]) * (
            1
            - PROBABILITIES["SIFT"]
            .loc[PROBABILITIES["SIFT"]["Score"] == sift, "Deleterious"]
            .values[0]
        )

    if polyphen >= TEST_MINIMUM_CUTOFFS["PolyPhen"]:
        score += (polyphen / TEST_MAXIMUM_CUTOFFS["PolyPhen"]) * (
            1
            - PROBABILITIES["PolyPhen"]
            .loc[PROBABILITIES["PolyPhen"]["Score"] == polyphen, "Normal"]
            .values[0]
        )
    else:
        score += (polyphen / TEST_MAXIMUM_CUTOFFS["PolyPhen"]) * (
            1
            - PROBABILITIES["PolyPhen"]
            .loc[PROBABILITIES["PolyPhen"]["Score"] == polyphen, "Deleterious"]
            .values[0]
        )

    # This acocunts for number of VEP teests used. Since we are only using 2, we use 2...
    score = score / 2

    if score >= 0.469:
        pred = "deleterious"
    elif score >= 0 and score < 0.469:
        pred = "neutral"

    return (score, pred)


def fishers_exact_test(
    input_row: dict, reference_population: str, comparison_population: list[str]
):
    """Runs a row-wise Fisher's Exact Test between the two listed populations
    Parameters:
    input (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data.
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """
    for fishers_exact_internal_key, fishers_exact_internal_dataset in input_row.items():
        fishers_exact_internal_dataset["OR"] = fishers_exact_internal_dataset["Count"][
            ["ID", "POS", "REF", "ALT"]
        ]
        fishers_exact_internal_dataset["P"] = fishers_exact_internal_dataset["Count"][
            ["ID", "POS", "REF", "ALT"]
        ]
        for population in comparison_population:
            for key2, direction in {
                "L": "less",
                "G": "greater",
                "T": "two-sided",
            }.items():
                o_label = f"{reference_population}_{key2}_{population}"
                p_label = f"{reference_population}_{key2}_{population}"
                fishers_exact_internal_dataset["OR"][o_label] = None
                fishers_exact_internal_dataset["P"][p_label] = None
                print(fishers_exact_internal_dataset["Count"].iterrows())
                for (
                    for_loop_var_index,
                    for_loop_var_row,
                ) in fishers_exact_internal_dataset["Count"].iterrows():
                    contingency = [
                        [
                            for_loop_var_row[f"{reference_population}_ac"],
                            for_loop_var_row[f"{reference_population}_tc"]
                            - for_loop_var_row[f"{reference_population}_ac"],
                        ],
                        [
                            for_loop_var_row[f"{population}_ac"],
                            for_loop_var_row[f"{population}_tc"]
                            - for_loop_var_row[f"{population}_ac"],
                        ],
                    ]
                    o_val, p_val = fisher_exact(contingency, alternative=direction)
                    fishers_exact_internal_dataset["P"].loc[
                        for_loop_var_index, p_label
                    ] = p_val
                    fishers_exact_internal_dataset["OR"].loc[
                        for_loop_var_index, o_label
                    ] = o_val
                fishers_exact_internal_dataset["P"][p_label] = multipletests(
                    fishers_exact_internal_dataset["P"][p_label], method="bonferroni"
                )[1]
        # columnsToDrop = list()
        # for pop in populations:
        #     columnsToDrop.append("{}_ac".format(pop))
        #     columnsToDrop.append("{}_tc".format(pop))
        # dataset.drop(columns=columnsToDrop ,inplace=True)


# %%
#  Import Data:
data = dict()
for gene in LOCATION_NAMES:
    data[gene] = read_vcf(join("..", "..", "results", "FINAL", f"ALL_{gene}.vcf.gz"))

# %%
# Sub-Divide data:
data_generator = dict()
for dataset in data.keys():
    data_generator[dataset] = chunk(data[dataset], 100)


# %%
# Compile and format request bodies:
data_to_send = dict()
# Iterate through each gene:
for dataset in LOCATION_NAMES:
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
        REQUESTING = True
        temp_list = list()
        while REQUESTING:
            r = requests.post(
                ENDPOINT,
                headers=HEADERS,
                data=json.dumps(chunk),
                params=generate_params(dataset_key),
            )
            if not r.ok:
                print(str(r.reason))
                time.sleep(2)
            else:
                REQUESTING = False
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
                transcripts_requested = TRANSCRIPTS.query(
                    f"gene_name == '{dataset_key}'"
                )["transcript_id"].tolist()
                for transcript in transcripts_requested:
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
                            s, p = condel_weighted_score(
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
                            s, p = condel_weighted_score(
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
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "Supplementary Table", cluster))
    for gene in LOCATION_NAMES:
        supplementary[gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_VEP.csv",
            ),
            sep="\t",
            index=False,
        )
        # supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)


# FREQUENCY CALCULATIONS


# Retrive VEP data:
supplementary = dict()

for cluster in CLUSTERS:
    supplementary[cluster] = dict()
    for gene in LOCATION_NAMES:
        supplementary[cluster][gene] = pd.read_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_VEP.csv",
            ),
            sep="\t",
        )[["ID", "POS", "REF", "ALT"]]


# Calculate frequencies
frequency_data = dict()
fishers_data = dict()


for cluster in CLUSTERS:
    frequency_data[cluster] = dict()
    fishers_data[cluster] = dict()
    for gene in LOCATION_NAMES:
        fishers_data[cluster][gene] = supplementary[cluster][gene][
            ["ID", "POS", "REF", "ALT"]
        ]
        for pop in CLUSTER_DATA[cluster].unique():
            supplementary[cluster][gene][pop] = 0
            fishers_data[cluster][gene][f"{pop}_ac"] = 0
            fishers_data[cluster][gene][f"{pop}_tc"] = 0
            frequency_data[cluster][gene] = pd.read_csv(
                join(
                    "..",
                    "..",
                    "results",
                    "FINAL",
                    cluster,
                    f"ALL_{gene}.{pop}.acount",
                ),
                delimiter="\t",
            ).rename(columns={"#CHROM": "CHROM"})
            for index, row in frequency_data[cluster][gene].iterrows():
                supplementary[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"], pop
                ] = freq(row["ALT_CTS"], row["OBS_CT"])
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    f"{pop}_ac",
                ] = row["ALT_CTS"]
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    f"{pop}_tc",
                ] = row["OBS_CT"]


# Save the resulting dataframe back to its excel file:
for cluster in CLUSTERS:
    for gene in LOCATION_NAMES:
        supplementary[cluster][gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_Freq.csv",
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
                f"{gene}_Count.csv",
            ),
            index=False,
            sep="\t",
        )


# FISHERS EXACT


# Load the Supplementary Table
supplementary = dict()

for cluster in CLUSTERS:
    supplementary[cluster] = dict()
    for gene in LOCATION_NAMES:
        supplementary[cluster][gene] = dict()
        supplementary[cluster][gene]["Count"] = pd.read_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_Count.csv",
            ),
            sep="\t",
        )
        supplementary[cluster][gene]["OR"] = dict()
        supplementary[cluster][gene]["P"] = dict()


# Run Fisher's Exact test
fishers_exact_test(supplementary["SUPER"], REFERENCE_POPULATION, COMPARISON_POPULATIONS)


# Save Fishers Data to CSV
for gene in LOCATION_NAMES:
    supplementary["SUPER"][gene]["P"].to_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            "SUPER",
            f"{gene}_FishersP.csv",
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
            f"{gene}_FishersOR.csv",
        ),
        sep="\t",
        index=False,
    )

# %%
