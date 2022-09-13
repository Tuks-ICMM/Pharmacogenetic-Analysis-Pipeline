#!/usr/bin/env python
"""
A Python script designed to run Fishers Exact Test with Bonferonni corrections and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from functools import reduce
from os.path import join
from common.common import directory_exists

import pandas as pd
from pandas import read_excel, merge

from common.common import save_or_append_to_excel

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
CLUSTERS.remove("SUB")

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# [SET] Statistical tests:
TEST_OUTPUTS = ["P", "OR"]
SHEET_NAMES = [
    "VEP",
    "Freq",
    "Count",
    "Fishers_two_sided",
]

# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        DATA[cluster][gene] = dict()
        for sheet in SHEET_NAMES:
            DATA[cluster][gene][sheet] = read_excel(
                join("..", "..", "results", "FINAL", f"{cluster}-{gene}.xlsx"),
                sheet_name=sheet,
            )


# %%
############ PERFORM JOIN ############
######################################
for cluster in CLUSTERS:
    for gene in GENES:
        ANALYSES = list()
        for sheet in SHEET_NAMES:
            ANALYSES.append(DATA[cluster][gene][sheet])
        DATA[cluster][gene]["ALL"] = reduce(
            lambda left, right: merge(
                left, right, on=["CHROM", "POS", "ID", "REF", "ALT"], how="inner"
            ),
            ANALYSES,
        )

# %%
############ SAVE TO EXCEL ############
#######################################

# [SAVE] VEP analysis to excel
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        save_or_append_to_excel(DATA[cluster][gene]["ALL"], cluster, gene, "ALL", "replace")
# %%
