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

from pandas import read_csv

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
    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################
    # [IMPORT] Data and initial data-templates with variant IDs and population columns:
    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
    )[MULTIINDEX]
    DATA.rename(columns={"#CHROM": "CHROM"})
    DATA.set_index(MULTIINDEX)
    _logger.info("Variant index data has been imported.")

    # [FOR] each dataset to be merged in
    for data in snakemake.input.analyses:
        # [IMPORT] dataset to be merged
        DATASET = read_csv(data, on_bad_lines="warn")
        DATASET.set_index(MULTIINDEX)
        _logger.info("Dataset successfully imported: '%s'", data) 

        DATA = DATA.merge(DATASET, how="left", left_index=True, right_index=True)
        _logger.info("Dataset successfully merged in: '%s'", DATA)

    _logger.info("Exporting data to output...")
    DATA.to_csv(snakemake.output.consolidated_data)
    _logger.info("Exported data to output...")

except Exception as E:
    _logger.error(E)
