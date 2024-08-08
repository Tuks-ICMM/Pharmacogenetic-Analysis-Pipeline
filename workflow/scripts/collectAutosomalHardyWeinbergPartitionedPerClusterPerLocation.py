#!/usr/bin/env python
"""
A Python script designed to run Frequency calculations and save the results to an Excel file.

[NEEDS] variant_count.py
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import logging

from pandas import read_csv, merge
from re import search

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

    _logger.info("VARIABLES: CLUSTER,LOCATION")
    _logger.info("CLUSTER: %s", snakemake.wildcards.cluster)
    _logger.info("LOCATION: %s", snakemake.wildcards.location)

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] Data and initial data-templates with variant IDs and population columns:
    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
        comment="#",
        sep="\t",
        names=MULTIINDEX + ["QUAL", "FILTER"],
        dtype={
            "CHROM": "Int8",
            "POS": "Int64",
        },
    )[MULTIINDEX]
    _logger.info("Variant index data has been imported.")

    ########################################################
    ############ COLLECT HARDY WEINBERG RESULTS ############
    ########################################################

    _logger.info(
        "Gathering Hardy-Weinberg reports generated. The following reports were found: %s",
        snakemake.input.hardy_weinberg,
    )
    for report in snakemake.input.hardy_weinberg:
        POPULATION = search("hardy_weinberg.([a-zA-Z]+).hardy", report).group(1)
        _logger.info("The %s population has been identified for analysis.", POPULATION)

        _logger.info("Importing report: %s", report)
        REPORT = read_csv(
            report,
            sep="\t",
            dtype={
                "#CHROM": "Int8",
                "POS": "Int64",
            },
        ).rename(
            columns={
                "#CHROM": "CHROM",
                "A1": f"{POPULATION}_A1",
                "HOM_A1_CT": f"{POPULATION}_HOM_A1_CT",
                "HET_A1_CT": f"{POPULATION}_HET_A1_CT",
                "TWO_AX_CT": f"{POPULATION}_TWO_AX_CT",
                "O(HET_A1)": f"{POPULATION}_O(HET_A1)",
                "E(HET_A1)": f"{POPULATION}_E(HET_A1)",
                "MIDP": f"{POPULATION}_MIDP",
            }
        )
        _logger.info("Imported report: %s", report)

        # _logger.info("Setting multi-index.")
        # REPORT.set_index(MULTIINDEX, inplace=True)
        # _logger.info("Set multi-index successfully.")

        DATA = merge(
            DATA,
            REPORT,
            how="left",
            on=MULTIINDEX,
            suffixes=("", "_clash"),
        )
        _logger.info(
            "Report has been filed and processed: %s",
            report,
        )
    DATA.to_csv(snakemake.output[0], index=False)
    _logger.info("Results written to file.")

except Exception as E:
    _logger.error(E)
