#!/usr/bin/env python
"""
A Python script designed to run Fishers Exact Test with Bonferonni corrections and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import join

import pandas as pd
from pandas import read_excel

from common.common import directory_exists, save_or_append_to_excel

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
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


# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        DATA[cluster][gene] = read_excel(
            join("..", "..", "results", "FINAL", f"{cluster}-{gene}.xlsx"),
            sheet_name="ALL",
        )

# %%
############ VEP AND ALLELE SHEET ############
##############################################
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        data_to_save = DATA[cluster][gene].query(
            " | ".join(
                [
                    f"(Diseases.notnull() & {population} >= 0.01)"
                    for population in ALL_POPULATIONS
                ]
            )
        )
        save_or_append_to_excel(data_to_save, cluster, gene, "VEP_ALLELE", "replace")
# %%
############ ALLELE AND FISHERS  SHEET ############
#####################################################
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        FISHERS_FILTER = " | ".join(
            [
                f"{REFERENCE_POPULATION}_P_{population} <= 0.05"
                for population in COMPARISON_POPULATIONS
            ]
        )
        data_to_save = DATA[cluster][gene].query(
            " | ".join(
                [
                    f"( {population} >= 0.01 & ({FISHERS_FILTER}))"
                    for population in ALL_POPULATIONS
                ]
            )
        )
        save_or_append_to_excel(
            data_to_save, cluster, gene, "ALLELE_FISHERS", "replace"
        )
# %%
############ VEP,ALLELE AND FISHERS SHEET ############
#####################################################
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        FISHERS_FILTER = " | ".join(
            [
                f"{REFERENCE_POPULATION}_P_{population} <= 0.05"
                for population in COMPARISON_POPULATIONS
            ]
        )
        data_to_save = DATA[cluster][gene].query(
            " | ".join(
                [
                    f"( Diseases.notnull() & {population} >= 0.01 & ({FISHERS_FILTER}))"
                    for population in ALL_POPULATIONS
                ]
            )
        )
        save_or_append_to_excel(
            data_to_save, cluster, gene, "VEP_ALLELE_FISHERS", "replace"
        )
# %%
