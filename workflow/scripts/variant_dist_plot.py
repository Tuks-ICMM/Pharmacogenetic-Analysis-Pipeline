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

    # [IMPORT] variant count data:
    VEP_RESULTS = read_csv(
        snakemake.input.vep_results,
        on_bad_lines="warn",
    )
    _logger.info(
        "Variant effect prediction results for subgroup %s | %s  has been imported.",
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
    )
    _logger.debug(VEP_RESULTS)

    VEP_RESULTS.set_index(MULTIINDEX, inplace=True)
    _logger.info("Multiindex has been successfully set for VEP results.")
    _logger.debug(VEP_RESULTS)

    POPULATION_SUBSET = read_csv(
        snakemake.input.freq,
        on_bad_lines="warn",
        sep=",",
    )
    _logger.info(
        "Population-level variant index for the %s subgroup %s | %s  has been imported.",
        snakemake.wildcards.population,
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
    )
    _logger.debug(POPULATION_SUBSET)

    POPULATION_SUBSET.set_index(MULTIINDEX, inplace=True)
    _logger.info("Multiindex has been successfully set for subgroup.")
    _logger.debug(POPULATION_SUBSET)

    POPULATION_SUBSET.query(f"{snakemake.wildcards.population} > 0", inplace=True)
    _logger.info(
        "Variants not observed in the %s population have been removed.",
        snakemake.wildcards.population,
    )
    _logger.debug(POPULATION_SUBSET)

    POPULATION_SUBSET = POPULATION_SUBSET.merge(
        VEP_RESULTS, left_index=True, right_index=True
    )
    _logger.info(
        "The population-level data has been annotated with the corresponding variant-effect results."
    )

    # EXPLODED_DATA = POPULATION_SUBSET.explode(["Variant Impact"])
    # EXPLODED_DATA.reset_index(inplace=True)
    # _logger.debug(
    #     "Variant impacts have been successfully collected and expanded to row-per-record format."
    # )
    # _logger.debug(POPULATION_SUBSET)
    POPULATION_SUBSET.reset_index(inplace=True)

    PLOT = px.histogram(
        POPULATION_SUBSET,
        x="POS",
        color="Impact",
        marginal="box",
        title=f"{snakemake.wildcards.location} | {snakemake.wildcards.cluster} -> {snakemake.wildcards.population} Variant distribution by impact",
        template="simple_white",
        color_discrete_map={
            "High": "red",
            "Moderate": "orange",
            "Low": "blue",
            "Modifier": "Purple",
        },
        category_orders={"Variant Impact": ["Low", "Moderate", "High", "Modifier"]},
        nbins=20,
    )
    _logger.debug("Plot has been successfully drafted.")

    PLOT.update_layout(margin=dict(l=0, r=0, b=5))
    _logger.debug("Plot margins have been configured.")

    PLOT.update_xaxes(tickformat="r")
    _logger.debug("Plot x-axis label formats have been configured.")

    PLOT.write_image(snakemake.output[0], scale=2)
    _logger.debug("Plot has been written to output.")

except Exception as E:
    _logger.error(E)
