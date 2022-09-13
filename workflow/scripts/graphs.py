#!/usr/bin/env python
"""
A Python script designed to generate an UpSet plot and save the resultq.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import json
import re
from itertools import cycle, islice
from os.path import join

import myvariant as mv
import numpy as np
import pandas as pd
from seaborn import set_style
from matplotlib.pyplot import savefig, suptitle, rcParams
from numpy import nan
from pandas import read_excel
from upsetplot import plot

set_style("whitegrid")
rcParams.update({"font.size": 7})

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


def ALL_POPULATIONS(cluster_level: str) -> list:
    return SAMPLE_DATA[cluster_level].unique().tolist()


# [SET] Statistical tests:
TEST_OUTPUTS = ["P", "OR"]
ALTERNATIVE_TESTS = [
    "less",
    "greater",
    "two-sided",
]


# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
# We want this to run per-cluster:
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    POPULATIONS = ALL_POPULATIONS(cluster)
    ALLELE_FILTERS = " | ".join([f"{population} >= 0.01" for population in POPULATIONS])
    # We want this to run per-gene:
    for gene in GENES:
        # [1] Import the data in
        DATA[cluster][gene] = read_excel(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}.xlsx",
            ),
            sheet_name="Freq",
        ).query(ALLELE_FILTERS)
        for population in POPULATIONS:
            DATA[cluster][gene][population] = DATA[cluster][gene][population].astype(
                bool
            )
        DATA[cluster][gene] = (
            DATA[cluster][gene].set_index(POPULATIONS).groupby(level=POPULATIONS).size()
        )


# %%
############ GENERATE GRAPH ############
########################################

# We want this to run per-cluster:
for cluster in CLUSTERS:
    # We want this to run per-gene:
    for gene in GENES:
        plot(
            DATA[cluster][gene],
            sort_by="cardinality",
            sort_categories_by="cardinality",
            show_percentages=True,
        )
        # pup.catplot(kind='')
        suptitle(
            f"{gene} alleles (>=1%) accros population intersections",
            fontweight="bold",
            fontsize="xx-large",
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
