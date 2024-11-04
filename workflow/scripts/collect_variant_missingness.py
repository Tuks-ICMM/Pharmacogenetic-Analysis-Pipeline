#!/usr/bin/env python
"""
A Python script designed to collect the count data from the Plink-2.0 reports. SAves to an Excel file.


[NEEDS] WORKFLOW_RUNTIME

[NEEDED_BY] fishers_exact_bonferonni_corrected.py
[NEEDED_BY] frequency_calculations.py

"""

############ IMPORT DEPENDANCIES ############
#############################################


import logging
from re import search
import traceback

from pandas import read_csv, merge

from scripts.common.common import MULTIINDEX

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


try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
        comment="#",
        sep="\t",
        dtype={
            "CHROM": "Int8",
            "POS": "Int64",
            "ID": "string",
            "REF": "string",
            "ALT": "string",
            "QUAL": "string",
            "FILTER": "string",
        },
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
    )
    _logger.info("Variant index data has been imported.")
    _logger.info(DATA)

    _logger.info(
        "Beginning population variant-missingness report files from available files found: [%s].",
        snakemake.input.missingness_reports,
    )

    for report in snakemake.input.missingness_reports:

        # Identify reported population.
        POPULATION = search(
            f"{snakemake.wildcards.cluster}/{snakemake.wildcards.location}/missingness_per_cluster/missingness.([a-zA-Z]+).vmiss.zst",
            report,
        ).group(1)
        _logger.info(
            "The %s population allele count report has been identified for analysis.",
            POPULATION,
        )

        REPORT = read_csv(
            report,
            header=0,
            on_bad_lines="warn",
            sep="\t",  # This imports the columns under the following names. Length equals number of columns
            dtype={
                "#CHROM": "Int8",
                "POS": "Int64",
                "ID": "string",
                "REF": "string",
                "ALT": "string",
                "PROVISIONAL_REF?": "string",
                "MISSING_CT": "Int64",
                "OBS_CT": "Int64",
                "F_MISS": "Float64",
            },
        ).rename(
            columns={
                "#CHROM": "CHROM",
                "PROVISIONAL_REF?": f"{POPULATION}_provisional",
                "MISSING_CT": f"{POPULATION}_obs_miss",
                "OBS_CT": f"{POPULATION}_obs_count",
                "F_MISS": f"{POPULATION}_f_miss",
            }
        )
        _logger.info("Imported report successfully:")
        _logger.info(REPORT)

        DATA = merge(DATA, REPORT, on=MULTIINDEX, how="left")
        
    DATA.drop(columns=["QUAL", "FILTER"], inplace=True)

    DATA.to_csv(snakemake.output.variant_missingness_report, index=False)
except Exception as E:
    _logger.error(E)
    _logger.error(traceback.format_exc())
    exit(1)
