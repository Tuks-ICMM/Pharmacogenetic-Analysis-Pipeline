#!/usr/bin/env python
"""
A Python script designed to generate a ternary plot to illustrate the Hardy-Weinberg results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import join

import numpy as np
import pandas as pd
import plotly.express as px
from matplotlib.pyplot import rcParams, savefig
from numpy import nan
from pandas import read_csv
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

DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        DATA[cluster][gene] = dict()
        for population in SAMPLE_DATA[cluster].unique():
            DATA[cluster][gene][population] = read_csv(
                join(
                    "..",
                    "..",
                    "results",
                    "FINAL",
                    cluster,
                    f"ALL_{gene}.{population}.hardy",
                ),
                sep="\t",
                dtype={
                    "#CHROM": str,
                    "ID": str,
                    "A1": str,
                    "AX": str,
                    "HOM_A1_CT": int,
                    "HET_A1_CT": int,
                    "TWO_AX_CT": int,
                    "O(HET_A1)": float,
                    "E(HET_A1)": float,
                    "MIDP": float,
                },
            )

# %%


for cluster in CLUSTERS:
    for gene in GENES:
        for population in SAMPLE_DATA[cluster].unique():
            DATA[cluster][gene][population]["HOM_ALT_FREQ"] = (
                DATA[cluster][gene][population]["HOM_A1_CT"]
                / (
                    DATA[cluster][gene][population]["HOM_A1_CT"]
                    + DATA[cluster][gene][population]["HET_A1_CT"]
                    + DATA[cluster][gene][population]["TWO_AX_CT"]
                )
                * 100
            )
            DATA[cluster][gene][population]["HET_FREQ"] = (
                DATA[cluster][gene][population]["HET_A1_CT"]
                / (
                    DATA[cluster][gene][population]["HOM_A1_CT"]
                    + DATA[cluster][gene][population]["HET_A1_CT"]
                    + DATA[cluster][gene][population]["TWO_AX_CT"]
                )
                * 100
            )
            DATA[cluster][gene][population]["HOM_REF_FREQ"] = (
                DATA[cluster][gene][population]["TWO_AX_CT"]
                / (
                    DATA[cluster][gene][population]["HOM_A1_CT"]
                    + DATA[cluster][gene][population]["HET_A1_CT"]
                    + DATA[cluster][gene][population]["TWO_AX_CT"]
                )
                * 100
            )
            DATA[cluster][gene][population]["HWE_P_VALUE_CUTOFF"] = np.where(
                DATA[cluster][gene][population]["MIDP"] <= 0.05,
                "Significant",
                "Not Significant",
            )

# %%


for cluster in CLUSTERS:
    for gene in GENES:
        for population in SAMPLE_DATA[cluster].unique():
            fig = px.scatter_ternary(
                DATA[cluster][gene][population],
                a="HET_FREQ",
                b="HOM_ALT_FREQ",
                c="HOM_REF_FREQ",
                color="HWE_P_VALUE_CUTOFF",
                labels={
                    "HWE_P_VALUE_CUTOFF": "Deviation from Hardy-Weinberg",
                    "HET_FREQ": "Heterozygous Freq.",
                    "HOM_ALT_FREQ": "Homozygous Alt. Freq.",
                    "HOM_REF_FREQ": "Homozygous Ref. Freq.",
                },
                title=f"{population} population {gene} gene Hardy-Weinburg (mid-p adjusted) results",
            )
            # savefig(
            #     join(
            #         "..",
            #         "..",
            #         "results",
            #         "FINAL",
            #         cluster,
            #         f"HWE_{gene}_{population}.jpeg",
            #     ),
            #     dpi=1200,
            # )
            fig.show()

# %%
