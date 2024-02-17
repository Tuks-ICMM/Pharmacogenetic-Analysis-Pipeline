#!/usr/bin/env python
"""
A Python script designed to run Fishers Exact Test with Bonferonni corrections and save the results to an Excel file.

[NEEDS] variant_count.py
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import join

import pandas as pd
from common.common import save_or_append_to_excel
from common.fishers_exact_test import fishers_exact_test
from numpy import nan
from pandas import Series, read_excel
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

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

# # [SET] Statistical test options needed for underlying statistical packages:
# TEST_OUTPUTS = ["P", "OR"]

# [SET] Types of Fishers Exact Tests to run.
ALTERNATIVE_TESTS = [
    "less",
    "greater",
    "two-sided",
]

# [SET] Multi-indexe columns
MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
for alternative_test_name in ALTERNATIVE_TESTS:
    DATA[alternative_test_name] = dict()
    for cluster in CLUSTERS:
        DATA[alternative_test_name][cluster] = dict()
        for gene in GENES:
            DATA[alternative_test_name][cluster][gene] = read_excel(
                join("..", "..", "results", "FINAL", f"{cluster}-{gene}.xlsx"),
                sheet_name="Count",
            )

            # for population in SAMPLE_DATA[cluster].unique():
            #     DATA[alternative_test_name][cluster][gene].rename(
            #         columns={
            #             f"{population}_tc": f"{alternative_test_name}_{population}_tc",
            #             f"{population}_ac": f"{alternative_test_name}_{population}_ac",
            #         },
            #         inplace=True,
            #     )
            DATA[alternative_test_name][cluster][gene].set_index(
                MULTIINDEX, inplace=True
            )


# %%
############ SHAPE DATA AND PERFORM ANALYSIS ############
#########################################################

# We want this to run per-alternative-statistical-test:
for alternative_test_name in ALTERNATIVE_TESTS:
    # We want to run this per-cluster:
    for cluster in CLUSTERS:
        # We want to run this per-gene:
        for gene in GENES:
            # We want to run this per-comarison-population:
            for population in COMPARISON_POPULATIONS:
                # Generating new columns for this slice:
                P_VALUE_COLUMN = (
                    f"{alternative_test_name}.{REFERENCE_POPULATION}_P_{population}"
                )
                ODDS_RATIO_COLUMN = (
                    f"{alternative_test_name}.{REFERENCE_POPULATION}_OR_{population}"
                )
                REFERENCE_POPULATION_AC_COLUMN = f"{REFERENCE_POPULATION}_ac"
                REFERENCE_POPULATION_TC_COLUMN = f"{REFERENCE_POPULATION}_tc"
                ALTERNATE_POPULATION_AC_COLUMN = f"{population}_ac"
                ALTERNATE_POPULATION_TC_COLUMN = f"{population}_tc"

                # Use Pandas magic to create contingency table
                DATA[alternative_test_name][cluster][gene]["contingency"] = DATA[
                    alternative_test_name
                ][cluster][gene].apply(
                    lambda row: [
                        [
                            row[REFERENCE_POPULATION_AC_COLUMN],
                            row[REFERENCE_POPULATION_TC_COLUMN]
                            - row[REFERENCE_POPULATION_AC_COLUMN],
                        ],
                        [
                            row[ALTERNATE_POPULATION_AC_COLUMN],
                            row[ALTERNATE_POPULATION_TC_COLUMN]
                            - row[ALTERNATE_POPULATION_AC_COLUMN],
                        ],
                    ],
                    axis=1,
                )

                # Here, we perform the Fishers Exact test per row and unpack the p-value and odds-ratio into seperate, labled columns
                (
                    DATA[alternative_test_name][cluster][gene][ODDS_RATIO_COLUMN],
                    DATA[alternative_test_name][cluster][gene][P_VALUE_COLUMN],
                ) = zip(
                    *DATA[alternative_test_name][cluster][gene].apply(
                        lambda row: tuple(
                            fisher_exact(
                                row["contingency"], alternative=alternative_test_name
                            )
                        ),
                        axis=1,
                    )
                )
# %%
############ DROP UNDEEDED COLUMNS ############
###############################################

# We want this to run per-alternative-statistical-test:
for alternative_test_name in ALTERNATIVE_TESTS:
    # We want to run this per-cluster:
    for cluster in CLUSTERS:
        COLS_TO_REMOVE = [
            f"{population}_{count_type}"
            for population in SAMPLE_DATA[cluster].unique()
            for count_type in ["tc", "ac"]
        ]
        COLS_TO_REMOVE.append("contingency")
        # We want to run this per-gene:
        for gene in GENES:
            try:
                DATA[alternative_test_name][cluster][gene].drop(
                    labels=COLS_TO_REMOVE, inplace=True, axis=1
                )
            except Exception as e:
                print(e)
# %%
############ SAVE FISHERS EXACT TEST DATA TO FILE ############
##############################################################
for alternative_test_name in ALTERNATIVE_TESTS:
    for cluster in CLUSTERS:
        for gene in GENES:
            save_or_append_to_excel(
                DATA[alternative_test_name][cluster][gene].reset_index(),
                cluster,
                gene,
                f"Fishers_{alternative_test_name.replace('-', '_')}",
                "replace",
            )

# %%
