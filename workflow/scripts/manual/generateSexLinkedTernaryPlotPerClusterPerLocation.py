#!/usr/bin/env python
"""
A Python script designed to generate a ternary plot for a given Plink-2 hardy-weinberg autosome report


[NEEDS] WORKFLOW_RUNTIME


"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################


import logging
from re import search

import plotly.graph_objects as go
from numpy import where
from pandas import read_csv
from plotly.subplots import make_subplots

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Fatima Barmania",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

# %%

LOGGER = logging.getLogger("variant_count.py")
LOGGER.setLevel(logging.DEBUG)

LOG_FILE = logging.FileHandler(snakemake.log[0])
LOG_FILE.setLevel(logging.DEBUG)
LOGGER.addHandler(LOG_FILE)

# %%

try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    LOGGER.debug("BEGIN")
    for HARDY_WEINBERG_REPORT in snakemake.input.hardy_weinberg_report:
        POPULATION = search(
            "hardy_weinberg.([a-zA-Z]+).hardy.x", HARDY_WEINBERG_REPORT
        ).group(1)

        # LOGGER.debug(f"The population has been identified: {POPULATION}")

        HWE_DATA = read_csv(
            HARDY_WEINBERG_REPORT,
            sep="\t",
        )
        LOGGER.debug("HWE report has been imported")

        # [CONVERT] genotype-counts into genotype-frequencies:
        HWE_DATA["FEMALE_HET_FREQ"] = HWE_DATA["FEMALE_HET_A1_CT"] / (
            HWE_DATA["FEMALE_HOM_A1_CT"]
            + HWE_DATA["FEMALE_HET_A1_CT"]
            + HWE_DATA["FEMALE_TWO_AX_CT"]
        )
        LOGGER.debug("The female heterozygous frequency column has been calculated")
        HWE_DATA["FEMALE_HOM_ALT_FREQ"] = HWE_DATA["FEMALE_HOM_A1_CT"] / (
            HWE_DATA["FEMALE_HOM_A1_CT"]
            + HWE_DATA["FEMALE_HET_A1_CT"]
            + HWE_DATA["FEMALE_TWO_AX_CT"]
        )
        LOGGER.debug(
            "The female homozygous alternate frequency column has been calculated"
        )
        HWE_DATA["FEMALE_HOM_REF_FREQ"] = HWE_DATA["FEMALE_TWO_AX_CT"] / (
            HWE_DATA["FEMALE_HOM_A1_CT"]
            + HWE_DATA["FEMALE_HET_A1_CT"]
            + HWE_DATA["FEMALE_TWO_AX_CT"]
        )
        LOGGER.debug(
            "The female homozygous reference frequency column has been calculated"
        )

        HWE_DATA["MALE_ALT_FREQ"] = HWE_DATA["MALE_A1_CT"] / (
            HWE_DATA["MALE_A1_CT"] + HWE_DATA["MALE_AX_CT"]
        )
        LOGGER.debug("The male alternate frequency column has been calculated")
        HWE_DATA["MALE_NON_ALT_FREQ"] = HWE_DATA["MALE_A1_CT"] / (
            HWE_DATA["MALE_A1_CT"] + HWE_DATA["MALE_AX_CT"]
        )
        LOGGER.debug("The male reference frequency column has been calculated")

        # [BUILD] graph-object using data
        TERNARY_PLOT = make_subplots(
            rows=1, cols=2, subplot_titles=["Diploid", "Haploid"]
        )
        LOGGER.debug("The graph figure has been created.")
        TERNARY_PLOT.update_layout(
            {
                "ternary": {
                    "sum": 100,
                    "aaxis": {
                        "tickangle": 0,
                        "title": {"text": "Homozygous Alt. Freq."},
                    },
                    "baxis": {
                        "tickangle": 45,
                        "title": {"text": "Homozygous Ref. Freq."},
                    },
                    "caxis": {
                        "tickangle": -45,
                        "title": {"text": "Heterozygous Freq."},
                    },
                    "domain": {"row": 1, "col": 1},
                },
                ""
                "title": {
                    "text": f"{snakemake.wildcards.cluster} | {POPULATION} | {snakemake.wildcards.location} Hardy-Weinberg Test (Sex-Linked)"
                },
                "width": 1400,
                "height": 700,
                "scattermode": "overlay",
            }
        )
        LOGGER.debug("The graph has been configured (Title, etc)")

        FEMALE = HWE_DATA.loc[
            HWE_DATA["FEMALE_HOM_A1_CT"]
            == 0 & HWE_DATA["FEMALE_HET_A1_CT"]
            == 0 & HWE_DATA["FEMALE_TWO_AX_CT"]
            == 0
        ]
        FEMALE_NONSIGNIFICANT = FEMALE.loc[FEMALE["MIDP"] <= 0.05]
        FEMALE_SIGNIFICANT = FEMALE.loc[FEMALE["MIDP"] > 0.05]
        LOGGER.debug("The female samples have been partitioned")

        MALE = HWE_DATA.loc[HWE_DATA["MALE_A1_CT"] == 0 & HWE_DATA["MALE_AX_CT"] == 0]
        MALE_NONSIGNIFICANT = MALE.loc[MALE["MIDP"] <= 0.05]
        MALE_SIGNIFICANT = MALE.loc[MALE["MIDP"] > 0.05]
        LOGGER.debug("The male samples have been partitioned")

        TERNARY_PLOT.add_trace(
            go.Scatterternary(
                {
                    "name": "Significant deviation",
                    "mode": "markers",
                    "a": FEMALE_NONSIGNIFICANT["HET_FREQ"].tolist(),
                    "b": FEMALE_NONSIGNIFICANT["HOM_ALT_FREQ"].tolist(),
                    "c": FEMALE_NONSIGNIFICANT["HOM_REF_FREQ"].tolist(),
                    "marker": {"symbol": "circle", "color": "grey", "opacity": 0.75},
                }
            ),
            row=1,
            col=1,
        )
        LOGGER.debug("Female variants with insignificant HWE deviation have been added")

        TERNARY_PLOT.add_trace(
            go.Scatter(
                {
                    "name": "Significant deviation",
                    "mode": "markers",
                    "a": MALE_NONSIGNIFICANT["MALE_ALT_FREQ"].tolist(),
                    "b": MALE_NONSIGNIFICANT["MALE_NON_ALT_FREQ"].tolist(),
                    "marker": {"symbol": "circle", "color": "grey", "opacity": 0.75},
                }
            ),
            row=1,
            col=2,
        )
        LOGGER.debug("Male variants with insignificant HWE deviation have been added")

        TERNARY_PLOT.add_trace(
            go.Scatterternary(
                {
                    "name": "Insignificant deviation",
                    "mode": "markers",
                    "a": FEMALE_SIGNIFICANT["HET_FREQ"].tolist(),
                    "b": FEMALE_SIGNIFICANT["HOM_ALT_FREQ"].tolist(),
                    "c": FEMALE_SIGNIFICANT["HOM_REF_FREQ"].tolist(),
                    "marker": {"symbol": "circle", "color": "red", "opacity": 0.75},
                }
            )
        )
        LOGGER.debug("Female variants with significant HWE deviation have been added")

        TERNARY_PLOT.add_trace(
            go.Scatter(
                {
                    "name": "Insignificant deviation",
                    "mode": "markers",
                    "a": MALE_SIGNIFICANT["MALE_ALT_FREQ"].tolist(),
                    "b": MALE_SIGNIFICANT["MALE_NON_ALT_FREQ"].tolist(),
                    "marker": {"symbol": "circle", "color": "red", "opacity": 0.75},
                }
            ),
            row=1,
            col=2,
        )
        LOGGER.debug("Male variants with significant HWE deviation have been added")

        TERNARY_PLOT.write_image(snakemake.output[0])
        LOGGER.debug("The graph has been written to file.")

# %%
except Exception as E:
    LOGGER.error(E)
