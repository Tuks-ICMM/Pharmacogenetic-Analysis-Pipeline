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
import traceback

from pandas import read_csv
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

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

############ CONFIGURE LOGGING ############
###########################################

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


############ SET CONSTANTS ############
#######################################

# [SET] Multi-indexe columns
MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]


# %%
try:
    ############ IMPORT DATA ############
    #####################################
    # [IMPORT] variant count data:
    DATA = read_csv(snakemake.input.allele_counts, on_bad_lines="warn", sep=",")
    _logger.info(
        "Variant count data for the %s subgroup has been imported successfully. Available keys: %s",
        snakemake.wildcards.location,
        ", ".join(DATA.keys()[:-1]) + f"{DATA.keys()[-1]}",
    )

    # [IMPORT] Sample annotations
    ANNOTATIONS = read_csv(snakemake.input.psam, sep="\t")
    _logger.info(
        "Sample annotations for the %s subgroup have been imported successfully.",
        snakemake.wildcards.location,
    )

    ############ SHAPE DATA AND PERFORM ANALYSIS ############
    #########################################################

    # [SET] Reference population
    REFERENCE_POPULATION_AC_COLUMN = f"{snakemake.params.reference_population}_ac"
    REFERENCE_POPULATION_TC_COLUMN = f"{snakemake.params.reference_population}_tc"
    _logger.debug(
        "Reference population %s has been set for analysis.",
        snakemake.params.reference_population,
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

    # We want to run this per-comarison-population:
    for population in COMPARISON_POPULATIONS:
        # [GENERATE] new columns for this population comparison:
        P_VALUE_COLUMN = f"{snakemake.params.reference_population}_P_{population}"
        ODDS_RATIO_COLUMN = f"{snakemake.params.reference_population}_OR_{population}"

        ALTERNATE_POPULATION_AC_COLUMN = f"{population}_ac"
        ALTERNATE_POPULATION_TC_COLUMN = f"{population}_tc"
        _logger.debug("The %s subgroup has been selected for analysis.", population)

        # [GENERATE] a contingency table in the form of an array (Pandas magic) that can be stored as-is in a single column
        DATA.apply(
            lambda row: _logger.debug(row[REFERENCE_POPULATION_AC_COLUMN]), axis=1
        )
        DATA["contingency"] = DATA.apply(
            lambda row: [
                [
                    row[REFERENCE_POPULATION_AC_COLUMN],
                    row[REFERENCE_POPULATION_TC_COLUMN]
                    - row[REFERENCE_POPULATION_AC_COLUMN],
                ],
                [
                    row[ALTERNATE_POPULATION_AC_COLUMN],
                    row[ALTERNATE_POPULATION_TC_COLUMN]
                    - row[ALTERNATE_POPULATION_AC_COLUMN],
                ],
            ],
            axis=1,
        )
        _logger.debug("Contingency table constructed")

        # [CALCULATE] Two-sided fishers exact test on each row unpack the p-value and odds-ratio into separate, labeled columns
        (
            DATA[ODDS_RATIO_COLUMN],
            DATA[P_VALUE_COLUMN],
        ) = zip(
            *DATA.apply(
                lambda row: tuple(
                    fisher_exact(row["contingency"], alternative="two-sided")
                ),
                axis=1,
            )
        )
        _logger.debug("Fishers_Exact calculation applied")

        # [ADJUST] P-values with Bonferoni multiple test corrections:
        (
            DATA[f"{population}_hypothesis_rejection"],
            DATA[f"{P_VALUE_COLUMN}_corrected"],
        ) = multipletests(DATA[P_VALUE_COLUMN], alpha=0.05, method="bonferroni")[:2]

        # [DELETE] Un-needed comparison population columns
        DATA.drop(
            labels=[ALTERNATE_POPULATION_AC_COLUMN, ALTERNATE_POPULATION_TC_COLUMN],
            inplace=True,
            axis=1,
        )
        _logger.debug(
            f"Redundant comparison population ({population}) table removed from memory"
        )

        # [DELETE] Un-needed contingency-table column
        DATA.drop(
            labels=["contingency"],
            inplace=True,
            axis=1,
        )
        _logger.debug("Temporary contingency table dropped")

    # [DELETE] Un-needed reference population column
    DATA.drop(
        labels=[REFERENCE_POPULATION_AC_COLUMN, REFERENCE_POPULATION_TC_COLUMN],
        inplace=True,
        axis=1,
    )
    _logger.debug(
        f"Redundant reference population allele counts dropped ({snakemake.params.reference_population})"
    )

    # %%
    ############ SAVE COUNT DATA TO FILE ############
    #################################################
    DATA.to_csv(snakemake.output[0], index=False)
    _logger.debug("Results saved to output")
# %%

except Exception as E:
    _logger.error(E)
    _logger.error(traceback.format_exc())
    exit(1)
