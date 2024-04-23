#!/usr/bin/env python
"""
A Python script designed to collect the count data from the Plink-2.0 reports. SAves to an Excel file.


[NEEDS] WORKFLOW_RUNTIME

[NEEDED_BY] fishers_exact_bonferonni_corrected.py
[NEEDED_BY] frequency_calculations.py

"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import logging
import sys
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

# %%


# def countTranscriptConsequences(apiResponse: dict) -> int:
#     if "transcript_consequences" in apiResponse.keys():
#         return len(apiResponse["transcript_consequences"])
#     else:
#         return 0


# def collectConsequenceTypes(apiResponse: dict) -> List[str]:
#     CONSEQUENCE_TYPES = list()
#     if "transcript_consequences" in apiResponse.keys():
#         for consequence in apiResponse["transcript_consequences"]:
#             CONSEQUENCE_TYPES.extend(consequence["consequence_terms"])
#     else:
#         pass
#     return CONSEQUENCE_TYPES


# def collectBiotypes(variant_effect_prediction: dict) -> list[str]:
#     BIOTYPES = list()
#     if "transcript_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["transcript_consequences"]:
#             BIOTYPES.append(consequence["biotype"].replace("_", " ").capitalize())
#     if "regulatory_feature_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["regulatory_feature_consequences"]:
#             BIOTYPES.append(consequence["biotype"])
#     return BIOTYPES


# def collectImpact(variant_effect_prediction: dict) -> list[str]:
#     IMPACT_PREDICTIONS = list()
#     if "transcript_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["transcript_consequences"]:
#             IMPACT_PREDICTIONS.append(
#                 consequence["impact"].replace("_", " ").capitalize()
#             )
#     if "regulatory_feature_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["regulatory_feature_consequences"]:
#             IMPACT_PREDICTIONS.append(
#                 consequence["impact"].replace("_", " ").capitalize()
#             )
#     return IMPACT_PREDICTIONS


# def collectFeatureTypes(variant_effect_prediction: dict) -> list[str]:
#     FEATURES = list()
#     if "transcript_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["transcript_consequences"]:
#             FEATURES.append("Transcript")
#     if "regulatory_feature_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["regulatory_feature_consequences"]:
#             FEATURES.append("Regulatory Feature")
#     return FEATURES


# %%

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
    _logger.debug(
        "Variant effect prediction results for subgroup %s | %s  has been imported.",
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
    )

    VEP_RESULTS.set_index(MULTIINDEX, inplace=True)
    _logger.debug("Multiindex has been successfully set.")
    _logger.debug(VEP_RESULTS)

    # [IMPORT] the frequency data:
    POPULATION_SUBSET = read_csv(
        snakemake.input.freq, on_bad_lines="warn", sep=",", header=0
    )
    _logger.debug(
        "Frequency results for subgroup %s -> %s -> %s  has been imported. Available keys are: [%s]",
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
        snakemake.wildcards.population,
        ", ".join(POPULATION_SUBSET.keys().to_list()),
    )

    POPULATION_SUBSET.query(f"{snakemake.wildcards.population} > 0", inplace=True)
    _logger.info(
        "Variants not observed in the %s population have been removed.",
        snakemake.wildcards.population,
    )
    _logger.debug(POPULATION_SUBSET)

    POPULATION_SUBSET.set_index(MULTIINDEX, inplace=True)
    _logger.info("Multiindex has been successfully set for subgroup.")
    _logger.debug(POPULATION_SUBSET)

    # [MERGE] using left-hand merge to combine and filter population-level VEP results for graphing:
    POPULATION_VEP_RESULTS = POPULATION_SUBSET.merge(
        VEP_RESULTS,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("_acount", None),
    )
    _logger.debug(
        "The Variant-Effect query results have been filtered and consolidated against the sub-population index."
    )
    _logger.debug(POPULATION_VEP_RESULTS)

    # %%

    # %%

    # DATA["transcript_count"] = DATA["result"].apply(lambda result: countTranscriptConsequences(result))

    # POPULATION_VEP_RESULTS["Variant Biotype"] = POPULATION_VEP_RESULTS["result"].apply(
    #     lambda result: collectBiotypes(result)
    # )
    # _logger.debug("Variant biotypes have been successfully collected.")

    # POPULATION_VEP_RESULTS["Variant Impact"] = POPULATION_VEP_RESULTS["result"].apply(
    #     lambda result: collectImpact(result)
    # )
    # _logger.debug("Variant impacts have been successfully collected.")

    # POPULATION_VEP_RESULTS["Feature Type"] = POPULATION_VEP_RESULTS["result"].apply(
    #     lambda result: collectFeatureTypes(result)
    # )
    # _logger.debug("Feature types have been successfully collected.")
    # _logger.debug(POPULATION_VEP_RESULTS)

    # EXPLODED_DATA = POPULATION_VEP_RESULTS.explode(
    #     ["Variant Biotype", "Variant Impact", "Feature Type"]
    # )
    # _logger.debug(
    #     "Variant biotypes, Variant imacts and Feature types have been successfully exploded into single rows for graphing."
    # )
    # _logger.debug(POPULATION_VEP_RESULTS)

    PLOT = px.sunburst(
        POPULATION_VEP_RESULTS.groupby(
            ["Feature_type", "Impact", "Biotype"],
            as_index=False,
        ).count(),
        path=["Feature_type", "Impact", "Biotype"],
        values="query",
        color="Impact",
        title=f"{snakemake.wildcards.cluster} | {snakemake.wildcards.location} | {snakemake.wildcards.population} Variant effect predictions",
        template="simple_white",
        color_discrete_map={
            "High": "red",
            "Moderate": "orange",
            "Low": "blue",
            "Modifier": "Purple",
        },
    )
    _logger.debug("Plot has been successfully drafted.")

    PLOT.update_layout(margin=dict(l=0, r=0, b=5))
    _logger.debug("Plot margins have been configured.")

    PLOT.write_image(snakemake.output[0], scale=2)
    _logger.debug("Plot has been written to output.")

# %%

except Exception as E:
    _logger.error(E)
    sys.exit()
