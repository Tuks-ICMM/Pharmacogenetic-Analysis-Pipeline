#!/usr/bin/env python
"""
A Python script designed to generate an UpSet plot and save the resultq.
"""
from os.path import join

import pandas as pd

# %%
############ IMPORT DEPENDANCIES ############
#############################################
import plotly.graph_objects as go
from matplotlib.pyplot import rcParams
from pandas import read_excel
from seaborn import set_style

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
        fig = go.Figure()
        fig.update_layout(
            title=f"Gnomad Frequency comarpsion per variant ({gene})",
            xaxis_title="Calculated Frequency",
            yaxis_title="Recorded Frequency",
            width=1200,
            height=600,
            autosize=False,
        )
        for population in SAMPLE_DATA[cluster].unique():
            POPULATION_LABEL = (
                f"colocated_variants.frequencies_gnomad.{population.lower()}"
            )
            if POPULATION_LABEL in DATA[cluster][gene]:
                fig.add_scatter(
                    x=DATA[cluster][gene][population],
                    y=DATA[cluster][gene][POPULATION_LABEL],
                    name=population,
                    mode="markers",
                )
        fig.show()

        # fig.write_image(
        #     join(
        #         "..",
        #         "..",
        #         "results",
        #         "FINAL",
        #         f"LINEAR_REGRESSION_{cluster}-{gene}_allele.jpeg",
        #     ),
        #     # dpi=1200,
        # )

# %%
############ JointPlot ############
# TODO: Investigate wide-to-long transformations
# TODO: Investigate JointPlots for regressions
###################################

# # We want this to run per-cluster:
# for cluster in CLUSTERS:
#     # We want this to run per-gene:
#     for gene in GENES:
#         jointplot(DATA[cluster][gene])
