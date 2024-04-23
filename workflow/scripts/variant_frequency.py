#!/usr/bin/env python
"""
A Python script designed to run Frequency calculations and save the results to an Excel file.

[NEEDS] variant_count.py
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import logging
from json import load
from os.path import join

from common.common import calculate_frequency, read_vcf
from pandas import read_csv

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
try:
    ############ IMPORT DAT2A & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] Data and initial data-templates with variant IDs and population columns:
    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
        comment="#",
        sep="\t",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
    )[MULTIINDEX]
    _logger.info("Variant index data has been imported.")
    # _logger.debug(DA2TA)

    ############ CALCULATE VARIANT FREQUENCIES ############
    #######################################################

    # [READ] Allele-count report
    ALLELE_COUNTS = read_csv(snakemake.input.allele_counts)
    _logger.info(
        "The %s | %s allele count report has been imported.",
        snakemake.wildcards.location,
        snakemake.wildcards.cluster,
    )
    _logger.info(ALLELE_COUNTS)

    # [BUILD] list of populations to calculate frequencies for
    POPULATIONS = [
        column.replace("_tc", "") for column in ALLELE_COUNTS if column.endswith("_tc")
    ]
    _logger.info(
        "The following populations have been identified for analysis: [%s]", POPULATIONS
    )

    # [BUILD] Allele-frequency report from copy of allele-count report
    ALLELE_FREQUENCIES = ALLELE_COUNTS.copy(deep=True)

    # [REPEAT] for each population
    for population in POPULATIONS:
        _logger.info(
            "Calculating frequencies for the %s population.",
            population,
        )
        # [CALCULATE] allele frequencies for a population column-wise
        ALLELE_FREQUENCIES[population] = ALLELE_FREQUENCIES.apply(
            lambda row: calculate_frequency(
                row[f"{population}_ac"], row[f"{population}_tc"]
            ),
            axis=1,
            result_type="reduce",
        )
        _logger.info(
            "The frequency of all variants found in the %s population have been calculated.",
            population,
        )
        # [DELETE] the allele-count columns for this population (memory-management)
        ALLELE_FREQUENCIES.drop(
            [f"{population}_ac", f"{population}_tc"], axis=1, inplace=True
        )
        _logger.info("The Allele-Count data has been removed form the output.")

    ############ SAVE TO OUTPUT ############
    #######################################
    ALLELE_FREQUENCIES.to_csv(snakemake.output[0], index=False)

except Exception as E:
    _logger.error(E)
