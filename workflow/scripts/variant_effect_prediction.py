#!/usr/bin/env python
"""
A Python script designed to run Variant Effect Prediction calculations
and API calls as well as save the results to file.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from json import dumps
from os.path import join
from time import sleep

import pandas as pd
from requests import post

from common.common import (
    chunk,
    directory_exists,
    generate_notation,
    generate_params,
    merge,
    read_vcf,
    save_or_append_to_excel,
)
from common.condel_score import condel_weighted_score

# [SET] Pandas Chained-Assignment off to prevent writing changes to transient DF copies:
pd.set_option("chained_assignment", None)


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
############ IMPORT DATA ############
#####################################

# [IMPORT] sample-level data:
SAMPLE_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = pd.read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = pd.read_csv(join("..", "..", "input", "transcripts.csv"))


# %%
############ SET CONSTANT VARIABLES ############
################################################

# [SET] cluster names:
CLUSTERS = set(SAMPLE_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# [SET] variant effect prediction metadata
ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


# %%
############ READ IN VCF FILE AS BASE ############
##################################################

RAW_DATA_TEMPLATE = dict()
for gene in GENES:
    RAW_DATA_TEMPLATE[gene] = read_vcf(
        join("..", "..", "results", "FINAL", f"ALL_{gene}.vcf.gz")
    )[["CHROM", "POS", "ID", "REF", "ALT"]]
    RAW_DATA_TEMPLATE[gene]["CHROM"] = RAW_DATA_TEMPLATE[gene]["CHROM"].str.extract(
        "chr([1-9]{1,2}|[XY])"
    )

# %%
############ SUB_DIVIDE DATA ############
#########################################
DATA_IN_CHUNKS = dict()
for genomic_location_name, dataset in RAW_DATA_TEMPLATE.items():
    DATA_IN_CHUNKS[genomic_location_name] = chunk(
        RAW_DATA_TEMPLATE[genomic_location_name], 100
    )


# %%
############ COMPILE E! ENSEMBL REQUEST PAYLOADS ############
#############################################################
# Compile and format request bodies:
API_PAYLOADS = dict()
# Iterate through each gene:
for dataset in GENES:
    API_PAYLOADS[dataset] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in DATA_IN_CHUNKS[dataset]:
        temp_list = list()

        # Iterate through each row in the chunk and add the HGVS notation to the list:
        for index, row in chunk.iterrows():
            temp_list.extend(
                generate_notation(
                    row,
                    GENE_DATA.loc[
                        GENE_DATA["location_name"] == dataset, "strand"
                    ].item(),
                )
            )
        API_PAYLOADS[dataset].append(dict(variants=temp_list))

# %%
############ PERFORM E! ENSEMBL API CALLS ############
######################################################
API_RESPONSE = dict()
for genomic_location_name, dataset in API_PAYLOADS.items():
    API_RESPONSE[genomic_location_name] = list()
    for index, chunk in enumerate(dataset):
        REQUESTING = True
        temp_list = list()
        while REQUESTING:
            r = post(
                ENDPOINT,
                headers=HEADERS,
                data=dumps(chunk),
                params=generate_params(genomic_location_name),
            )
            if not r.ok:
                print(str(r.reason))
                sleep(2)
            else:
                decoded = r.json()
                API_RESPONSE[genomic_location_name] = (
                    API_RESPONSE[genomic_location_name] + decoded
                )
                REQUESTING = False


# %%
############ COMPILE E! ENSEMBL API RESPONSES ############
##########################################################
# Abstract response parsing to functions
##########################################################
DATA = dict()
for genomic_location_name, genomic_location_responses in API_RESPONSE.items():
    DATA[genomic_location_name] = RAW_DATA_TEMPLATE[genomic_location_name][
        ["CHROM", "ID", "POS", "REF", "ALT"]
    ]
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
        DATA[genomic_location_name][column] = None

    # for _, genomic_location in API_RESPONSE.items():
    for variant in genomic_location_responses:
        # row = supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] ==
        # int(variant['start'])]
        DATA[genomic_location_name].loc[
            DATA[genomic_location_name]["POS"] == int(variant["start"]),
            "Start Coordinates",
        ] = variant["start"]
        DATA[genomic_location_name].loc[
            DATA[genomic_location_name]["POS"] == int(variant["start"]), "input"
        ] = variant["input"]

        co_variants = list()
        if "colocated_variants" in variant:
            DATA[genomic_location_name].loc[
                DATA[genomic_location_name]["POS"] == int(variant["start"]),
                "Co-Located Variant",
            ] = True
            for colocated_variant in variant["colocated_variants"]:
                co_variants.append(colocated_variant["id"])
        else:
            DATA[genomic_location_name].loc[
                DATA[genomic_location_name]["POS"] == int(variant["start"]),
                "Co-Located Variant",
            ] = False
            co_variants.append("-")
        DATA[genomic_location_name].loc[
            DATA[genomic_location_name]["POS"] == int(variant["start"]),
            "Existing Variation",
        ] = "| ".join(co_variants)

        if "transcript_consequences" in variant:
            transcripts_requested = TRANSCRIPT_DATA.query(
                f"gene_name == '{genomic_location_name}'"
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

                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "SIFT_score",
                    ] = (
                        consequence["sift_score"]
                        if ("sift_score" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "SIFT_pred",
                    ] = (
                        consequence["sift_prediction"]
                        if ("sift_prediction" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Polyphen_score",
                    ] = (
                        str(consequence["polyphen_score"])
                        if ("polyphen_score" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Polyphen_pred",
                    ] = (
                        consequence["polyphen_prediction"]
                        if ("polyphen_prediction" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Diseases",
                    ] = merge(phenotype)
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Consequence",
                    ] = merge(consequence["consequence_terms"])
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Transcript ID",
                    ] = consequence["transcript_id"]
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Biotype",
                    ] = (
                        consequence["biotype"] if ("biotype" in consequence) else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "CADD_PHRED",
                    ] = (
                        consequence["cadd_phred"]
                        if ("cadd_phred" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Transcript Strand",
                    ] = (
                        consequence["strand"] if ("strand" in consequence) else None
                    )
                    # row['LoFtool'] = consequence['loftool'] if 'loftool' in consequence
                    # else None
                    if "sift_score" in consequence and "polyphen_score" in consequence:
                        s, p = condel_weighted_score(
                            int(consequence["sift_score"]),
                            consequence["polyphen_score"],
                        )
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL",
                        ] = s
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL_pred",
                        ] = p
                    else:
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL",
                        ] = None
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL_pred",
                        ] = None
                    break

                else:
                    consequence = variant["transcript_consequences"][0]
                    phenotype = set()

                    if "phenotypes" in consequence:
                        for instance in consequence["phenotypes"]:
                            phenotype.add(instance["phenotype"])

                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "SIFT_score",
                    ] = (
                        consequence["sift_score"]
                        if ("sift_score" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "SIFT_pred",
                    ] = (
                        consequence["sift_prediction"]
                        if ("sift_prediction" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Polyphen_score",
                    ] = (
                        str(consequence["polyphen_score"])
                        if ("polyphen_score" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Polyphen_pred",
                    ] = (
                        consequence["polyphen_prediction"]
                        if ("polyphen_prediction" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Diseases",
                    ] = merge(phenotype)
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Consequence",
                    ] = merge(consequence["consequence_terms"])
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Transcript ID",
                    ] = consequence["transcript_id"]
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Biotype",
                    ] = (
                        consequence["biotype"] if ("biotype" in consequence) else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "CADD_PHRED",
                    ] = (
                        consequence["cadd_phred"]
                        if ("cadd_phred" in consequence)
                        else None
                    )
                    DATA[genomic_location_name].loc[
                        DATA[genomic_location_name]["POS"] == int(variant["start"]),
                        "Transcript Strand",
                    ] = (
                        consequence["strand"] if ("strand" in consequence) else None
                    )
                    # row['LoFtool'] = consequence['loftool'] if 'loftool' in
                    # consequence else None
                    if "sift_score" in consequence and "polyphen_score" in consequence:
                        s, p = condel_weighted_score(
                            int(consequence["sift_score"]),
                            consequence["polyphen_score"],
                        )
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL",
                        ] = s
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL_pred",
                        ] = p
                    else:
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL",
                        ] = None
                        DATA[genomic_location_name].loc[
                            DATA[genomic_location_name]["POS"] == int(variant["start"]),
                            "CONDEL_pred",
                        ] = None
                    break

# %%
############ SAVE TO EXCEL ############
#######################################

# [SAVE] VEP analysis to excel
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        save_or_append_to_excel(DATA[gene], cluster, gene, "VEP", "replace")
# %%
