#!/usr/bin/env python
"""
A Python script designed to Combine all sheets in an Excel file into an '_all' sheet.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from functools import reduce
from os.path import join
from typing import List

import pandas as pd
from common.common import directory_exists, save_or_append_to_excel
from pandas import merge

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
############ IMPORT DATA & SET CONSTANTS ############
#####################################################

# [IMPORT] sample-level data:
SAMPLE_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = pd.read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = pd.read_csv(join("..", "..", "input", "transcripts.csv"))

# [SET] cluster names:
CLUSTERS = set(SAMPLE_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")
CLUSTERS.remove("sex")
CLUSTERS.remove("sub-population")

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]


# [SET] API payload categories for results:
VARIANT_TYPE_CATEGORIES = [
    "colocated_variants",
    "transcript_consequences",
    "regulatory_feature_consequences",
]


def removeFromList(
    sheet_names: List[int | str], names_to_remove: List[str]
) -> List[int | str]:
    """A function to remove unintended sheets from the excel before processing"""
    RESULT = sheet_names
    for sheet_name in names_to_remove:
        try:
            RESULT.remove(sheet_name)
        except Exception as e:
            print(e)
    return RESULT


# [SET] Statistical tests:
TEST_OUTPUTS = ["P", "OR"]
# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        EXCEL_FILE = pd.ExcelFile(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}.xlsx",
            )
        )
        DATA[cluster][gene] = dict()
        for sheet in EXCEL_FILE.sheet_names:
            DATA[cluster][gene][sheet] = EXCEL_FILE.parse(sheet)
            if "ID" in DATA[cluster][gene][sheet]:
                DATA[cluster][gene][sheet].drop("ID", inplace=True, axis=1)

# %%
############ PERFORM JOIN (NON VEP) ############
######################################
for cluster in CLUSTERS:
    print(f"START CLUSTER: {cluster}")
    for gene in GENES:
        print(f"START GENE: {gene}")
        ANALYSES = list()
        for sheet in removeFromList(list(DATA[cluster][gene].keys()), ["_all"]):
            # removeFromList(
            #     list(DATA[cluster][gene].keys())), VARIANT_TYPE_CATEGORIES + ["_all", "ID"]
            # ):
            print(f"START SHEET: {sheet}")
            ANALYSES.append(DATA[cluster][gene][sheet])
        DATA[cluster][gene]["_all"] = reduce(
            lambda left, right: merge(
                left, right, on=["CHROM", "POS", "REF", "ALT"], how="inner"
            ),
            ANALYSES,
        )

# %%
############ PERFORM JOIN (VEP) ############
############################################
# for cluster in CLUSTERS:
#     for gene in GENES:
#         ANALYSES = list()
#         for sheet in removeFromList(
#             list(DATA[cluster][gene].keys()),
#             removeFromList(list(DATA[cluster][gene].keys()), VARIANT_TYPE_CATEGORIES),
#         ):
#             ANALYSES.append(DATA[cluster][gene][sheet].drop("ID", inplace=True))
#         DATA[cluster][gene]["_all"] = pd.merge(
#             DATA[cluster][gene]["_all"],
#             reduce(
#                 lambda left, right: merge(
#                     left, right, on=["CHROM", "POS", "REF", "ALT"], how="inner"
#                 ),
#                 ANALYSES,
#             ),
#             on=["CHROM", "POS", "REF", "ALT"],
#             how="inner",
#         )
# %%
############ SAVE TO EXCEL ############
#######################################

# [SAVE] VEP analysis to excel
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        save_or_append_to_excel(
            DATA[cluster][gene]["_all"], cluster, gene, "_all", "replace"
        )
# %%
