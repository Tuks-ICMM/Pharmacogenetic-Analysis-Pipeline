#!/usr/bin/env python
"""A Python script designed to merge multiple CSV files into an excel file.
"""

# Import Dependancies
import glob
import json
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
with open(join("config", "config.json")) as f:
    config = json.load(f)
genes = config["locations"]
clusters = config["cluster"]["clusters"]
populations = ["AFR", "AMR", "EUR", "EAS", "SAS"]
tests = ["VEP", "Freq", "Count", "FishersP", "FishersOR"]


# Import data to merge:
data = dict()
for file in glob.glob(join("results", "Supplementary Table", "*", "*_*.csv")):
    if platform.system() == "Windows":
        path = file.replace("\\", "/")
    else:
        path = file

    g = re.search("results/Supplementary Table/(.*)/(.*)_(.*).csv", path)
    cluster = g.group(1)
    gene = g.group(2)
    test = g.group(3)

    if not cluster in data:
        data[cluster] = dict()
        data[cluster][gene] = dict()
        data[cluster][gene][test] = dict()
    elif not gene in data[cluster]:
        data[cluster][gene] = dict()
        data[cluster][gene][test] = dict()
    elif not test in data[cluster][gene]:
        data[cluster][gene][test] = dict()

    data[cluster][gene][test] = pd.read_csv(
        join(
            "results",
            "Supplementary Table",
            g.group(1),
            "{gene}_{test}.csv".format(gene=gene, test=test),
        ),
        delimiter="\t",
    )


# Import data to merge:
# data = dict()
# for cluster in clusters:
#     data[cluster] = dict()
#     for gene in genes:
#         data[cluster][gene] = dict()
#         for test in tests:
# if cluster == 'SUB':
#     if (test != 'FishersP' and test != 'FishersOR'):
#         data[cluster][gene][test] = pd.read_csv("../final/Supplementary Table/{cluster}/{gene}_{test}.csv".format(cluster=cluster, gene=gene, test=test), delimiter='\t')
#     else:
#         pass
# else:
#     data[cluster][gene][test] = pd.read_csv("../final/Supplementary Table/{cluster}/{gene}_{test}.csv".format(cluster=cluster, gene=gene, test=test), delimiter='\t')


# Compile a VEP + Freq sheet for analysis sake:
for cluster in data.keys():
    for gene in data[cluster].keys():
        tests = list(data[cluster][gene].keys())
        data[cluster][gene]["ALL"] = data[cluster][gene][tests[0]][
            ["ID", "POS", "REF", "ALT"]
        ]
        try:
            tests.remove("ALL")
        except:
            pass
        for test in tests:
            if test == "FishersP" or test == "FishersOR":
                suffix = "_" + re.search("^Fishers([A-Z]{1,2})$", test).group(1)
                print("SUFFIX:", suffix)
                data[cluster][gene]["ALL"] = pd.merge(
                    data[cluster][gene]["ALL"],
                    data[cluster][gene][test],
                    suffixes=[None, suffix],
                    on=["ID", "POS", "REF", "ALT"],
                )
            else:
                data[cluster][gene]["ALL"] = pd.merge(
                    data[cluster][gene]["ALL"], data[cluster][gene][test]
                )


# Save data to Excel
for cluster in clusters:
    for gene in genes:
        with pd.ExcelWriter(
            join(
                "results",
                "{cluster}-{gene}.xlsx".format(cluster=cluster, gene=gene),
            ),
            engine="xlsxwriter",
        ) as writer:
            for test in data[cluster][gene]:
                df = data[cluster][gene][test]
                df.to_excel(writer, sheet_name=test, index=False, header=False)
                worksheet = writer.sheets["{}".format(test)]
                (max_row, max_col) = df.shape
                options = {
                    "columns": [{"header": column} for column in df.columns],
                    "name": test,
                    "first_column": True,
                }
                worksheet.add_table(0, 0, max_row - 1, max_col - 1, options)
