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
from matplotlib.pyplot import rcParams, savefig, suptitle
from numpy import nan
from pandas import read_excel
from seaborn import set_style
from upsetplot import UpSet

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
CLUSTERS.remove("sub-population")
CLUSTERS.remove("sex")

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
            sheet_name="_all",
        )
        for population in POPULATIONS:
            DATA[cluster][gene].set_index(
                DATA[cluster][gene][population] >= 0.01, append=True, inplace=True
            )
        DATA[cluster][gene] = DATA[cluster][gene].droplevel(None, axis=0)
        # DATA[cluster][gene] = (
        #     DATA[cluster][gene].set_index(POPULATIONS).groupby(level=POPULATIONS).size()
        # )


# %%
############ GENERATE Freq & VEP GRAPH ############
###################################################

# We want this to run per-cluster:
for cluster in CLUSTERS:
    # We want this to run per-gene:
    for gene in GENES:
        upset = UpSet(
            DATA[cluster][gene],
            subset_size="count",
            sort_by="cardinality",
            show_counts=True,
            sort_categories_by="cardinality",
            show_percentages=True,
            intersection_plot_elements=0,
            totals_plot_elements=4,
            min_degree=0.1,
        )
        # upset.add_stacked_bars(by="variant_class_y", title="Variant Class")
        upset.add_stacked_bars(
            by="most_severe_consequence_x", title="Most Severe Consequence"
        )
        upset.add_stacked_bars(by="transcript_consequences.biotype", title="Biotype")
        upset.add_stacked_bars(
            by="transcript_consequences.impact", title="Transcript Impact"
        )
        upset.add_stacked_bars(
            by="transcript_consequences.consequence_terms", title="Consequence Terms"
        )

        # [HIGHLIGHT] population specific partitions
        for population in POPULATIONS:
            ANTI_LIST = list(POPULATIONS)
            ANTI_LIST.remove(population)
            upset.style_subsets(present=population, absent=ANTI_LIST, facecolor="red")
            upset.style_subsets(present=POPULATIONS, facecolor="gray")
        suptitle(
            f"{cluster} {gene} alleles (>=1%) across population intersections",
            fontweight="bold",
            fontsize="xx-large",
        )
        upset.plot()
        savefig(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"UPSETPLOT_{cluster}-{gene}_allele.jpeg",
            ),
            dpi=1200,
        )
# %%
