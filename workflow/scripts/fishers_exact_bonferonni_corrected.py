#!/usr/bin/env python
"""
A Python script designed to run Fishers Exact Test with Bonferonni corrections and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import exists, join

import pandas as pd
from numpy import nan
from pandas import read_csv, read_excel

from common.common import read_vcf, save_or_append_to_excel
from common.fishers_exact_test import fishers_exact_test

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
ALTERNATIVE_TESTS = [
    "less",
    "greater",
    "two-sided",
]

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
                odds_ratio_label = f"{REFERENCE_POPULATION}_OR_{population}"
                p_value_label = f"{REFERENCE_POPULATION}_P_{population}"

                # Create the column
                DATA[alternative_test_name][cluster][gene][odds_ratio_label] = nan
                DATA[alternative_test_name][cluster][gene][p_value_label] = nan

                # We want this to repeat per-row (Hence why we need to correct later):
                for row in DATA[alternative_test_name][cluster][gene].to_dict(
                    "records"
                ):
                    results = fishers_exact_test(
                        row, REFERENCE_POPULATION, population, alternative_test_name
                    )

                    # Set the Odds Ratio:
                    DATA[alternative_test_name][cluster][gene].loc[
                        (
                            DATA[alternative_test_name][cluster][gene]["CHROM"]
                            == row["CHROM"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["ID"]
                            == row["ID"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["REF"]
                            == row["REF"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["ALT"]
                            == row["ALT"]
                        ),
                        odds_ratio_label,
                    ] = results["OR"]

                    # Set the P-value:
                    DATA[alternative_test_name][cluster][gene].loc[
                        (
                            DATA[alternative_test_name][cluster][gene]["CHROM"]
                            == row["CHROM"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["ID"]
                            == row["ID"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["REF"]
                            == row["REF"]
                        )
                        & (
                            DATA[alternative_test_name][cluster][gene]["ALT"]
                            == row["ALT"]
                        ),
                        p_value_label,
                    ] = results["P"]

                # DATA[cluster][gene][population] = DATA[cluster][gene].apply(lambda row: fishers_exact_test(row[REFERENCE_POPULATION], row[population]), result_type='expand')
                # fishers_exact_test(, REFERENCE_POPULATION, COMPARISON_POPULATIONS)

# %%
############ DROP UNDEEDED COLUMNS ############
###############################################

# We want this to run per-alternative-statistical-test:
for alternative_test_name in ALTERNATIVE_TESTS:
    # We want to run this per-cluster:
    for cluster in CLUSTERS:
        # We want to run this per-gene:
        for gene in GENES:
            try:
                DATA[alternative_test_name][cluster][gene].drop(columns=[f"{population}_{count_type}" for population in ALL_POPULATIONS for count_type in ['tc', 'ac']], inplace=True)
            except:
                pass
# %%
############ SAVE COUNT DATA TO FILE ############
#################################################
for alternative_test_name in ALTERNATIVE_TESTS:
    for gene in GENES:
        save_or_append_to_excel(
            DATA[alternative_test_name][cluster][gene],
            cluster,
            gene,
            f"Fishers_{alternative_test_name.replace('-', '_')}",
            "replace",
        )

# %%
