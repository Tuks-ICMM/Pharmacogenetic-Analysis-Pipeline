#!/usr/bin/env python
"""
A Python script designed to generate a ternary plot for a given Plink-2 hardy-weinberg autosome report


[NEEDS] WORKFLOW_RUNTIME


"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################


import logging

import plotly.express as px
from pandas import read_csv, merge

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


# def calcAverageCaddPhredScores(variant_effect_prediction: dict) -> str:
#     SCORES = list()
#     if "transcript_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["transcript_consequences"]:
#             SCORES.append(consequence["cadd_phred"])
#     if "regulatory_feature_consequences" in variant_effect_prediction:
#         for consequence in variant_effect_prediction["regulatory_feature_consequences"]:
#             SCORES.append(consequence["cadd_phred"])
#     return mean(SCORES)


# def collectMostSevereConsequence(variant_effect_prediction: dict) -> str:
#     if "most_severe_consequence" in variant_effect_prediction:
#         return (
#             variant_effect_prediction["most_severe_consequence"]
#             .replace("_", " ")
#             .capitalize()
#         )


# %%

try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns:
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    _logger.debug("BEGIN")

    # [IMPORT] the frequency data:
    POPULATION_SUBSET = read_csv(
        snakemake.input.freq_report, on_bad_lines="warn", sep=","
    )
    POPULATION_SUBSET.set_index(MULTIINDEX, inplace=True)

    CONSOLIDATED_DATA = read_csv(
        snakemake.input.consolidated_data,
        on_bad_lines="warn",
    )
    _logger.debug(
        "Frequency results for subgroup %s -> %s -> %s  has been imported. Available keys are: [%s]",
        snakemake.wildcards.cluster,
        snakemake.wildcards.location,
        snakemake.wildcards.population,
        ", ".join(CONSOLIDATED_DATA.keys().to_list()),
    )

    CONSOLIDATED_DATA.set_index(MULTIINDEX, inplace=True)
    _logger.info("Multi-index has been set.")
    _logger.debug(CONSOLIDATED_DATA)

    # HWE_DATA = read_csv(
    #     snakemake.input.hardy_weinberg_report[0],
    #     sep="\t",
    # )
    # _logger.info("HWE report has been imported")
    # _logger.debug(HWE_DATA)

    # HWE_DATA.rename(columns={"#CHROM": "CHROM"}, inplace=True)
    # _logger.info("Column names have been cleaned up.")

    # HWE_DATA.set_index(MULTIINDEX, inplace=True)
    # _logger.info("Multiindex has been set.")
    # _logger.debug(HWE_DATA)

    POPULATION_SUBSET = merge(
        POPULATION_SUBSET,
        CONSOLIDATED_DATA,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("", "_consolidated"),
    )
    _logger.info(
        "Consolidated data/results have been merged against register of variants for the %s population.",
        snakemake.wildcards.population,
    )

    # [IMPORT] variant count data:
    VEP_RESULTS = read_csv(
        snakemake.input.vep_results,
        on_bad_lines="warn",
    )
    _logger.info(
        "Variant effect prediction results for subgroup %s | %s  has been imported.",
        snakemake.wildcards.location,
        snakemake.wildcards.cluster,
    )
    _logger.debug(VEP_RESULTS)

    VEP_RESULTS.set_index(MULTIINDEX, inplace=True)
    _logger.info("The Multi-index has been created successfully.")
    _logger.debug(VEP_RESULTS)

    # VEP_RESULTS["result"] = VEP_RESULTS["result"].apply(lambda result: loads(result))
    # _logger.debug(
    #     "The VEP outputs have been successfully converted from JSON to in-memory python dict."
    # )

    # VEP_RESULTS["Av. CADD PHRED score"] = VEP_RESULTS["result"].apply(
    #     lambda result: calcAverageCaddPhredScores(result)
    # )
    # _logger.debug(
    #     "The VEP outputs have been successfully converted from JSON to in-memory python dict."
    # )

    # VEP_RESULTS["Most Severe Consequence"] = VEP_RESULTS["result"].apply(
    #     lambda result: collectMostSevereConsequence(result)
    # )
    # _logger.debug(
    #     "The VEP outputs have been successfully converted from JSON to in-memory python dict."
    # )

    POPULATION_SUBSET = merge(
        POPULATION_SUBSET,
        VEP_RESULTS,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("", "_vep"),
    )
    _logger.info(
        "VEP report has been consolidated against register of variants for the %s population.",
        snakemake.wildcards.population,
    )

    # # [IMPORT] variant count data:
    # FISHERS_RESULTS = read_csv(
    #     snakemake.input.fishers,
    #     on_bad_lines="warn",
    #     sep=",",
    # )
    # _logger.info(
    #     "Variant effect prediction results for subgroup %s | %s  has been imported.",
    #     snakemake.wildcards.location,
    #     snakemake.wildcards.cluster,
    # )
    # _logger.debug(FISHERS_RESULTS)

    # FISHERS_RESULTS.set_index(MULTIINDEX, inplace=True)
    # _logger.info("The Multi-index has been created successfully.")
    # _logger.debug(FISHERS_RESULTS)

    # POPULATION_SUBSET = POPULATION_SUBSET.merge(
    #     FISHERS_RESULTS,
    #     how="left",
    #     left_index=True,
    #     right_index=True,
    #     suffixes=("", "_fishers"),
    # )
    # _logger.info(
    #     "VEP report has been consolidated against register of variants for the %s population.",
    #     snakemake.wildcards.population,
    # )

    # [CONVERT] genotype-counts into genotype-frequencies:
    for population in snakemake.params.populations:
        pass
    POPULATION_SUBSET["Hetrozygous Freq."] = (
        POPULATION_SUBSET["HET_A1_CT"]
        / (
            POPULATION_SUBSET["HOM_A1_CT"]
            + POPULATION_SUBSET["HET_A1_CT"]
            + POPULATION_SUBSET["TWO_AX_CT"]
        )
    ) * 100
    _logger.debug("The heterozygous frequency column has been calculated")
    POPULATION_SUBSET["Homozygous Alt Freq."] = (
        POPULATION_SUBSET["HOM_A1_CT"]
        / (
            POPULATION_SUBSET["HOM_A1_CT"]
            + POPULATION_SUBSET["HET_A1_CT"]
            + POPULATION_SUBSET["TWO_AX_CT"]
        )
    ) * 100
    _logger.debug("The homozygous alternate frequency column has been calculated")
    POPULATION_SUBSET["Homozygous Ref Freq."] = (
        POPULATION_SUBSET["TWO_AX_CT"]
        / (
            POPULATION_SUBSET["HOM_A1_CT"]
            + POPULATION_SUBSET["HET_A1_CT"]
            + POPULATION_SUBSET["TWO_AX_CT"]
        )
    ) * 100
    _logger.debug("The homozygous reference frequency column has been calculated")

    # [BUILD] graph-object using data
    # TERNARY_PLOT = go.Figure()
    POPULATION_SUBSET.loc[POPULATION_SUBSET["MIDP"] > 0.05, "HWE Test"] = (
        "Non-Significant"
    )
    POPULATION_SUBSET.loc[POPULATION_SUBSET["MIDP"] <= 0.05, "HWE Test"] = (
        POPULATION_SUBSET.loc[POPULATION_SUBSET["MIDP"] <= 0.05, "Consequence_type"]
    )

    # POPULATION_SUBSET.loc[
    #     POPULATION_SUBSET[
    #         f"{snakemake.params.reference_population}_P_{snakemake.wildcards.population}"
    #     ]
    #     > 0.05,
    #     "Non-significant",
    # ] = False
    # POPULATION_SUBSET.loc[
    #     POPULATION_SUBSET[
    #         f"{snakemake.params.reference_population}_P_{snakemake.wildcards.population}"
    #     ]
    #     <= 0.05,
    #     "Significantly different frequency",
    # ] = True
    TERNARY_PLOT = px.scatter_ternary(
        POPULATION_SUBSET,
        a="Hetrozygous Freq.",
        b="Homozygous Alt Freq.",
        c="Homozygous Ref Freq.",
        text="query",
        color="HWE Test",
        # symbol="Significantly different frequency",
        # size=f"{snakemake.params.reference_population}_OR_{snakemake.wildcards.population}",
        # size="CADD_PHRED",
        template="simple_white",
        width=1280,
        height=960,
        color_discrete_map={"Non-Significant": "grey"},
        title=f"{snakemake.wildcards.location}  | {snakemake.wildcards.cluster} -> {population} Hardy-Weinberg Test",
    )
    _logger.debug("The graph figure has been created.")
    # TERNARY_PLOT.update_layout(
    #     {
    #         "ternary": {
    #             "sum": 100,
    #             "aaxis": {
    #                 "tickangle": 0,
    #                 "ticksuffix": " %",
    #                 "titlefont": {"size": 20},
    #                 "tickfont": {"size": 15},
    #                 "title": {"text": "Homozygous Alt. Freq."},
    #             },
    #             "baxis": {
    #                 "tickangle": 45,
    #                 "ticksuffix": " %",
    #                 "titlefont": {"size": 20},
    #                 "tickfont": {"size": 15},
    #                 "title": {"text": "Homozygous Ref. Freq."},
    #             },
    #             "caxis": {
    #                 "tickangle": -45,
    #                 "ticksuffix": " %",
    #                 "titlefont": {"size": 20},
    #                 "tickfont": {"size": 15},
    #                 "title": {"text": "Heterozygous Freq."},
    #             },
    #         },
    #         "title": {
    #             "text": f"{snakemake.wildcards.location}  | {snakemake.wildcards.cluster} -> {POPULATION} Hardy-Weinberg Test",
    #             "xanchor": "center",
    #             "x": 0.5,
    #             "yanchor": "top",
    #             "y": 0.95,
    #             "font": {"size": 25},
    #         },
    #         "margin": {"t": 150, "r": 80, "b": 80, "l": 80},
    #         "showlegend": True,
    #         "width": 1280,
    #         "height": 960,
    #         "scattermode": "overlay",
    #     },
    #     template="simple_white",
    # )
    _logger.debug("The graph has been configured (Title, etc)")

    NONSIGNIFICANT = POPULATION_SUBSET.loc[POPULATION_SUBSET["MIDP"] < 0.05]
    # TERNARY_PLOT.add_trace(
    #     go.Scatterternary(
    #         {
    #             "name": "Significant deviation",
    #             "mode": "markers",
    #             "a": NONSIGNIFICANT["Hetrozygous Freq."].tolist(),
    #             "b": NONSIGNIFICANT["Homozygous Alt Freq."].tolist(),
    #             "c": NONSIGNIFICANT["Homozygous Ref Freq."].tolist(),
    #             "marker": {"symbol": "circle", "color": "red", "opacity": 0.75},
    #         }
    #     )
    # )
    _logger.debug("Variants with insignificant HWE deviation have been added")

    SIGNIFICANT = POPULATION_SUBSET.loc[POPULATION_SUBSET["MIDP"] >= 0.05]
    # TERNARY_PLOT.add_trace(
    #     go.Scatterternary(
    #         {
    #             "name": "Non-significant deviation",
    #             "mode": "markers",
    #             "a": SIGNIFICANT["Hetrozygous Freq."].tolist(),
    #             "b": SIGNIFICANT["Homozygous Alt Freq."].tolist(),
    #             "c": SIGNIFICANT["Homozygous Ref Freq."].tolist(),
    #             "marker": {"symbol": "circle", "color": "grey", "opacity": 0.75},
    #             "text": [
    #                 "Heterozygous Freq.",  # A-axis name
    #                 "Homozygous Alt. Freq.",  # B-axis name
    #                 "Homozygous Ref. Freq.",  # C-axis name
    #             ],
    #         }
    #     )
    # )
    _logger.debug("Variants with significant HWE deviation have been added")
    TERNARY_PLOT.write_image(snakemake.output[0])
    _logger.debug("The graph has been written to file.")

# %%
except Exception as E:
    _logger.error(E)
