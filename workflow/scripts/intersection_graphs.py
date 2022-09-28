#!/usr/bin/env python
"""
A Python script designed to run Fishers Exact Test with Bonferonni corrections and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import join

from matplotlib.pyplot import suptitle, savefig
from numpy import nan
from pandas import read_csv, read_excel
from upsetplot import plot

from common.common import save_or_append_to_excel
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
SAMPLE_DATA = read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = read_csv(join("..", "..", "input", "transcripts.csv"))

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
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        DATA[cluster][gene] = read_excel(
            join("..", "..", "results", "FINAL", f"{cluster}-{gene}.xlsx"),
            sheet_name="Freq",
        )

        DATA[cluster][gene] = DATA[cluster][gene].query(
            " | ".join([f"{population} >= 0.01" for population in ALL_POPULATIONS])
        )

        for population in ALL_POPULATIONS:
            DATA[cluster][gene][population] = DATA[cluster][gene][population].astype(
                bool
            )

        DATA[cluster][gene] = DATA[cluster][gene].set_index(ALL_POPULATIONS)

# %%
############ GENERATE GRAPHS ############
#########################################

for cluster in CLUSTERS:
    for gene in GENES:
        plot(
            DATA[cluster][gene],
            sort_by="cardinality",
            sort_categories_by="cardinality",
            show_percentages=True,
            subset_size="count",
        )
        suptitle(
            "{} alleles (>=1%) accros population intersections".format(gene),
            fontweight="semibold",
            fontsize="x-large",
        )
        savefig(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}_UpSetPlot.jpeg",
            ),
            dpi=1200,
        )

# %%
