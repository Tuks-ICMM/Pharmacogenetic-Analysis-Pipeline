#!/usr/bin/env python
"""
A Python script designed to generate a dist plot of all variant positions for a given location.


[NEEDS] WORKFLOW_RUNTIME

[NEEDED_BY] fishers_exact_bonferonni_corrected.py
[NEEDED_BY] frequency_calculations.py

"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################


import logging
from json import loads
from re import search
from typing import List

import plotly.express as px
from pandas import merge, read_csv

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


_logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=snakemake.log[0],
    encoding="utf-8",
    filemode="w",
    level=logging.DEBUG,
    format="[%(asctime)s] %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
_logger.info("---LOGGING INTERFACE STARTED---")


def collectImpact(variant_effect_prediction: dict) -> list[str]:
    IMPACT_PREDICTIONS = list()
    if "transcript_consequences" in variant_effect_prediction:
        for consequence in variant_effect_prediction["transcript_consequences"]:
            IMPACT_PREDICTIONS.append(
                consequence["impact"].replace("_", " ").capitalize()
            )
    if "regulatory_feature_consequences" in variant_effect_prediction:
        for consequence in variant_effect_prediction["regulatory_feature_consequences"]:
            IMPACT_PREDICTIONS.append(
                consequence["impact"].replace("_", " ").capitalize()
            )
    return IMPACT_PREDICTIONS


def collectMostSevereConsequence(variant_effect_prediction: dict) -> str:
    if "most_severe_consequence" in variant_effect_prediction:
        return (
            variant_effect_prediction["most_severe_consequence"]
            .replace("_", " ")
            .capitalize()
        )


try:
    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns:
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] Sample annotations
    ANNOTATIONS = read_csv(snakemake.input.psam, sep="\t")
    _logger.info(
        "Sample annotations for the %s subgroup have been imported successfully.",
        snakemake.wildcards.location,
    )

    # [CALCULATE] the list of comparison populations:
    COMPARISON_POPULATIONS: list = (
        ANNOTATIONS[snakemake.wildcards.cluster].unique().tolist()
    )
    COMPARISON_POPULATIONS.remove(snakemake.params.reference_population)
    _logger.debug(
        "Remaining populations set as comparison populations: [%s]",
        ", ".join(COMPARISON_POPULATIONS),
    )

    # [IMPORT] variant count data:
    FISHERS_RESULTS = read_csv(
        snakemake.input.fishers,
        on_bad_lines="warn",
        sep=",",
    )
    _logger.info(
        "Fishers-Exact results for subgroup %s | %s  has been imported.",
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
    )
    _logger.debug(FISHERS_RESULTS)

    FISHERS_RESULTS.set_index(MULTIINDEX, inplace=True)
    _logger.info("Multiindex has been successfully set for VEP results.")
    _logger.debug(FISHERS_RESULTS)

    # We want to run this per-comarison-population:
    FISHERS_RESULTS.query(
        " | ".join(
            [
                f"{snakemake.params.reference_population}_P_{population} <= 0.5"
                for population in COMPARISON_POPULATIONS
            ]
        ),
        inplace=True,
    )

    FISHERS_RESULTS = FISHERS_RESULTS.melt(
        value_vars=[
            f"{snakemake.params.reference_population}_P_{population}"
            for population in COMPARISON_POPULATIONS
        ]
        + [
            f"{snakemake.params.reference_population}_OR_{population}"
            for population in COMPARISON_POPULATIONS
        ],
        var_name="Comparisons",
        value_name="Score",
    )

    _logger.info(
        "Melting of Fishers Exact results is complete. Remaining keys are: [%s]",
        ", ".join(FISHERS_RESULTS.keys().tolist()),
    )
    _logger.debug(FISHERS_RESULTS)

    FISHERS_RESULTS.query(
        "Comparisons.str.contains('_OR_') & Score.notnull()", inplace=True
    )
    _logger.debug(FISHERS_RESULTS)

    PLOT = px.violin(
        FISHERS_RESULTS,
        x="Comparisons",
        y="Score",
        title="Significant Fishers-Exact test Odds-Ratio comparison",
        template="simple_white",
        box=True,
        # color="Variant Impact",
        # marginal="box",
        # title=f"{snakemake.wildcards.location} | {snakemake.wildcards.cluster} -> {snakemake.wildcards.population} Variant distribution by impact",
        # template="simple_white",
        # color_discrete_map={
        #     "High": "red",
        #     "Moderate": "orange",
        #     "Low": "blue",
        #     "Modifier": "Purple",
        # },
        # category_orders={"Variant Impact": ["Low", "Moderate", "High", "Modifier"]},
        # nbins=20,
    )
    _logger.debug("Plot has been successfully drafted.")

    PLOT.add_hline(
        y=1,
        line_dash="dot",
        annotation_text="Odds Ratio of 1",
        annotation_position="top right",
    )

    PLOT.update_layout(margin=dict(l=0, r=0, b=5))
    _logger.debug("Plot margins have been configured.")

    PLOT.update_xaxes(tickformat="r")
    _logger.debug("Plot x-axis label formats have been configured.")

    PLOT.write_image(snakemake.output[0], scale=2)
    _logger.debug("Plot has been written to output.")

except Exception as E:
    _logger.error(E)
