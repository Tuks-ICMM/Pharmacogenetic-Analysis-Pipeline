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
############ PERFORM VARIANT EFFECT PREDICTION ############
###########################################################

data = dict()
for gene in GENES:
    data[gene] = read_vcf(join("..", "..", "results", "FINAL", f"ALL_{gene}.vcf.gz"))

# %%
############ SUB_DIVIDE DATA ############
#########################################
data_generator = dict()
for dataset_key, dataset in data.items():
    data_generator[dataset_key] = chunk(data[dataset_key], 100)


# %%
############ COMPILE E! ENSEMBL REQUEST PAYLOADS ############
#############################################################
# Compile and format request bodies:
data_to_send = dict()
# Iterate through each gene:
for dataset in GENES:
    data_to_send[dataset] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in data_generator[dataset]:
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
        data_to_send[dataset].append(dict(variants=temp_list))

# %%
############ PERFORM E! ENSEMBL API CALLS ############
######################################################
data_received = dict()
for dataset_key, dataset in data_to_send.items():
    data_received[dataset_key] = list()
    for index, chunk in enumerate(dataset):
        REQUESTING = True
        temp_list = list()
        while REQUESTING:
            r = post(
                ENDPOINT,
                headers=HEADERS,
                data=dumps(chunk),
                params=generate_params(dataset_key),
            )
            if not r.ok:
                print(str(r.reason))
                sleep(2)
            else:
                REQUESTING = False
                decoded = r.json()
                data_received[dataset_key] = data_received[dataset_key] + decoded


# %%
############ COMPILE E! ENSEMBL API RESPONSES ############
##########################################################
# ToDo: Abstract response parsing to functions
##########################################################
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
            # row = supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] ==
            # int(variant['start'])]
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
                transcripts_requested = TRANSCRIPT_DATA.query(
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
                        # row['LoFtool'] = consequence['loftool'] if 'loftool' in consequence
                        # else None
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
                        # row['LoFtool'] = consequence['loftool'] if 'loftool' in
                        # consequence else None
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
############ SAVE TO EXCEL ############
#######################################

# [SAVE] VEP analysis to excel
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        supplementary[gene].to_excel(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}.xlsx",
            ),
            index=False,
            sheet_name="Variant-EFfect-Prediction",
        )
# START REDUNDANT DELETE
# supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)
# %%

# FREQUENCY CALCULATIONS

# %%
# Retrive VEP data:
# supplementary = dict()

# for cluster in CLUSTERS:
#     supplementary[cluster] = dict()
#     for gene in GENES:
#         supplementary[cluster][gene] = pd.read_csv(
#             join(
#                 "..",
#                 "..",
#                 "results",
#                 "Supplementary Table",
#                 f"{cluster}-{gene}.csv",
#             ),
#             sep="\t",
#         )[["ID", "POS", "REF", "ALT"]]

# END REDUNDANT DELETE
# %%


# %%
