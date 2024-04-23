#!/usr/bin/env python
"""
A Python script designed to collect the count data from the Plink-2.0 reports. SAves to an Excel file.


[NEEDS] WORKFLOW_RUNTIME

[NEEDED_BY] fishers_exact_bonferonni_corrected.py
[NEEDED_BY] frequency_calculations.py

"""

# pylint: disable=logging-fstring-interpolation
# %%
############ IMPORT DEPENDANCIES ############
#############################################


import logging
from json import loads
from sys import exit

from matplotlib import cm
from matplotlib.pyplot import savefig, suptitle
from pandas import read_csv
from upsetplot import UpSet

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


def collectMostSevereConsequence(variant_effect_prediction: dict) -> str:
    if "most_severe_consequence" in variant_effect_prediction:
        return (
            variant_effect_prediction["most_severe_consequence"]
            .replace("_", " ")
            .capitalize()
        )


# %%
try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] variant count data:
    DATA = read_csv(
        snakemake.input.freq,
        sep=",",
    )
    _logger.debug(
        f"Frequency data has been imported successfully for location {snakemake.wildcards.location}, available keys are: {DATA.keys()}"
    )

    DATA.rename(
        columns={
            population: f"{population}_freq"
            for population in snakemake.params.populations
        },
        inplace=True,
    )
    _logger.debug(
        f"Columns have been successfully renamed. Available indexes are now: {DATA.keys()}"
    )

    DATA.set_index(MULTIINDEX, inplace=True)
    _logger.info("The Multi-index has been created successfully.")
    _logger.debug(DATA)

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

    # [MERGE] in variant-effect-prediction results onto partition dataset:
    DATA = DATA.merge(
        VEP_RESULTS,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("", "_vep"),
    )
    _logger.info("Multi-index-based merge successful.")
    _logger.debug(DATA)

    # HWE_DATA = read_csv(
    #     snakemake.input.hardy_weinberg_report[0],
    #     sep="\t",
    # )
    # _logger.info("HWE report has been imported")
    # _logger.debug(HWE_DATA)

    # HWE_DATA.rename(columns={"#CHROM": "CHROM"}, inplace=True)
    # _logger.info("Column names have been cleaned up.")
    # HWE_DATA.loc[HWE_DATA["MIDP"] > 0.05, "Fishers significance"] = "Not significant"
    # HWE_DATA.loc[HWE_DATA["MIDP"] <= 0.05, "Fishers significance"] = "Significant"

    # HWE_DATA.set_index(MULTIINDEX, inplace=True)
    # _logger.info("Multiindex has been set.")
    # _logger.debug(HWE_DATA)

    # DATA = DATA.merge(
    #     HWE_DATA,
    #     how="left",
    #     left_index=True,
    #     right_index=True,
    #     suffixes=("", "_hwe"),
    # )
    # _logger.info("Multi-index-based merge successful.")
    # _logger.debug(DATA)

    for population in snakemake.params.populations:
        DATA.loc[(DATA[f"{population}_freq"] >= 0.01), population] = True
        DATA.loc[(DATA[f"{population}_freq"] < 0.01), population] = False
        _logger.debug(
            f"Conversion to boolean type completed for location: {population}"
        )

        # DATA.drop(columns=[f"{population}_freq"], inplace=True)
        # _logger.debug(f"Old Frequency column for {population} has been dropped.")

    _logger.info(
        f"Boolean conversion has been completed. Available keys include: {DATA.keys()}"
    )
    _logger.debug(DATA)

    DATA.query(
        " | ".join(
            [f"{population} >= 0.01" for population in snakemake.params.populations]
        ),
        inplace=True,
    )
    _logger.info(
        "Mutations that occurred at frequencies of less than 1% (i.e. rare variants) in any population group have been removed."
    )
    _logger.debug(DATA)

    _logger.debug(
        f"Attempting boolean multiindexing with available keys: {DATA.keys()} of length {len(DATA.keys())}"
    )
    _logger.debug(DATA)
    DATA.set_index(snakemake.params.populations, inplace=True)
    _logger.debug(
        f"Multi-index has been set successfully for populations: {snakemake.params.populations}."
    )
    _logger.debug(DATA)
    _logger.debug(DATA.groupby(by="Consequence_type").count())

    fig = UpSet(
        DATA,
        intersection_plot_elements=0,
        sort_by="cardinality",
        sort_categories_by="cardinality",
        show_percentages=True,
        subset_size="count",
        min_subset_size=DATA.shape[0]
        * 1
        / 100,  # TODO: Figure this out... at some point...
    )
    _logger.debug("UpPSetPlot has been successfully plotted.")

    fig.add_stacked_bars(
        by="Consequence_type",
        colors="tab20",
        title="Count by most severe consequence",
        elements=10,
    )
    _logger.info("Stacked bars for variant impact have been added onto the plot.")

    # fig.add_stacked_bars(
    #     by="Fishers significance",
    #     colors="tab20",
    #     title="Count by most severe consequence",
    #     elements=10,
    # )
    # _logger.info("Stacked bars for Fishers significance have been added onto the plot.")

    # fig.add_catplot(
    #     kind="violin",
    #     by="Variant Consequence",
    #     colors="tab20",
    #     title="Count by most severe consequence",
    #     elements=10,
    # )
    # _logger.info("Stacked bars for variant impact have been added onto the plot.")

    for population in snakemake.params.populations:
        fig.style_subsets(
            present=population,
            absent=[pop for pop in snakemake.params.populations if pop != population],
            edgecolor="red",
            linewidth=2,
        )
        _logger.info(
            "The %s partition has been outlined in red for reference purposes.",
            population,
        )

    fig.style_subsets(
        present=snakemake.params.populations, edgecolor="grey", linewidth=2
    )
    _logger.info(
        "The universal partition has been outlined in grey for reference purposes."
    )

    fig.plot()
    suptitle(
        f"{snakemake.wildcards.cluster} | {snakemake.wildcards.location} alleles (>=1%) across population intersections"
    )
    _logger.debug("Figure title has been configured successfully.")

    savefig(snakemake.output[0], dpi=100, bbox_inches="tight")
    _logger.debug("Plot has been saved successfully.")
# %%

except Exception as E:
    _logger.error(E)
    exit(E)
