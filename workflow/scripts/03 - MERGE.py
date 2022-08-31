#!/usr/bin/env python
"""
A Python script designed to merge multiple CSV files into an excel file.

This files sole purpose is to merge all available analysis results (.csv format)
into a single excel, using the names on the results files to determine where and
how the results should end up in the final Excel.
"""
# %%
# Import Dependancies
import glob
import platform
import re
from os.path import join

import pandas as pd

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Sarah Turner",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"


# Declare Constants and Functions:
LOCATIONS = pd.read_csv(join("..", "..", "input", "locations.csv"))
LOCATION_NAMES = LOCATIONS["location_name"].unique().tolist()

CLUSTER_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))
CLUSTERS = set(CLUSTER_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")

populations = ["AFR", "AMR", "EUR", "EAS", "SAS"]
tests = ["VEP", "Freq", "Count", "FishersP", "FishersOR"]


# Import data to merge, using iteration to import
# and compile each tests CSV results into a dict:
ANALYSIS_DATA = dict()
for file in glob.glob(
    join("..", "..", "results", "Supplementary Table", "*", "*_*.csv")
):
    if platform.system() == "Windows":
        path = file.replace("\\", "/")
    else:
        path = file

    # Determine the nature of the test based on the filename:
    g = re.search("results/Supplementary Table/(.*)/(.*)_(.*).csv", path)
    cluster = g.group(1)
    gene = g.group(2)
    test = g.group(3)

    # Format the placeholder dict according to each imported test:
    if not cluster in ANALYSIS_DATA:
        ANALYSIS_DATA[cluster] = dict()
        ANALYSIS_DATA[cluster][gene] = dict()
        ANALYSIS_DATA[cluster][gene][test] = dict()
    elif not gene in ANALYSIS_DATA[cluster]:
        ANALYSIS_DATA[cluster][gene] = dict()
        ANALYSIS_DATA[cluster][gene][test] = dict()
    elif not test in ANALYSIS_DATA[cluster][gene]:
        ANALYSIS_DATA[cluster][gene][test] = dict()

    # Actually import the data and save it to the placeholder variable:
    ANALYSIS_DATA[cluster][gene][test] = pd.read_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            g.group(1),
            f"{gene}_{test}.csv",
        ),
        delimiter="\t",
    )


# %%
# Compile a VEP + Freq sheet for analysis sake:
for data_key, cluster in ANALYSIS_DATA.items():
    for gene_key, gene in cluster.items():
        tests = list(gene.keys())
        ANALYSIS_DATA[data_key][gene_key]["ALL"] = ANALYSIS_DATA[data_key][gene_key][
            tests[0]
        ][["ID", "POS", "REF", "ALT"]]
        try:
            tests.remove("ALL")
        except Exception as e:  # pylint: disable=broad-except
            print("Error writing files: \n", e)
        for test in tests:
            if test == "FishersP" or test == "FishersOR":
                suffix = "_" + re.search("^Fishers([A-Z]{1,2})$", test).group(1)
                ANALYSIS_DATA[data_key][gene_key]["ALL"] = pd.merge(
                    ANALYSIS_DATA[data_key][gene_key]["ALL"],
                    ANALYSIS_DATA[data_key][gene_key][test],
                    suffixes=[None, suffix],
                    on=["ID", "POS", "REF", "ALT"],
                )
            else:
                ANALYSIS_DATA[data_key][gene_key]["ALL"] = pd.merge(
                    ANALYSIS_DATA[data_key][gene_key]["ALL"],
                    ANALYSIS_DATA[data_key][gene_key][test],
                )

# %%
# Save data to Excel
for cluster in CLUSTERS:
    for gene in LOCATION_NAMES:
        with pd.ExcelWriter(  # pylint: disable=abstract-class-instantiated
            join(
                "..",
                "..",
                "results",
                f"{cluster}-{gene}.xlsx",
            ),
            engine="xlsxwriter",
        ) as writer:
            for test in list(ANALYSIS_DATA[cluster][gene].keys().tolist()):
                data_to_write = ANALYSIS_DATA[cluster][gene][test]
                data_to_write.to_excel(
                    writer, sheet_name=test, index=False, header=False
                )
                worksheet = writer.sheets[f"{test}"]
                (row_count, col_count) = data_to_write.shape
                options = {
                    "columns": [{"header": column} for column in data_to_write.columns],
                    "name": test,
                    "first_column": True,
                }
                worksheet.add_table(0, 0, row_count - 1, col_count - 1, options)

# %%
