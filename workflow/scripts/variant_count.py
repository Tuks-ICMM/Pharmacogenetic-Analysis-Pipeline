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
from re import search

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
try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] index describing the variants to include. This will form the backbone of the output file and set variants to be worked with.
    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
        comment="#",
        sep="\t",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
    )[MULTIINDEX]
    _logger.info("Variant index data has been imported.")
    _logger.debug(DATA)
    # %%
    ############ MERGE RAW COUNTS ############
    ##########################################
    # [REPEAT] for each allele-count report produced per-population in this cluster

    _logger.info(
        "Beginning population count report files from available files found: [%s].",
        snakemake.input.allele_counts,
    )
    for report in snakemake.input.allele_counts:
        # [EXTRACT] the population code from the report name
        POPULATION = search("allele_count.([a-zA-Z]+).acount", report).group(1)
        _logger.info("The %s population has been identified for analysis.", POPULATION)

        # [IMPORT] allele-count report
        ALLELE_COUNTS = read_csv(
            report,
            on_bad_lines="warn",
            sep="\t",  # This imports the columns under the following names. Length equals number of columns
        )
        _logger.info(
            "The %s populations allele count reports have been imported.", POPULATION
        )
        ALLELE_COUNTS.rename(
            columns={
                "#CHROM": "CHROM",
                "ALT_CTS": f"{POPULATION}_ac",
                "OBS_CT": f"{POPULATION}_tc",
            },
            inplace=True,
        )
        ALLELE_COUNTS.drop(columns=["REF_CT"], inplace=True)

        # [MERGE(left)] allele-count report into base `Dataframe`
        DATA = merge(
            DATA,
            ALLELE_COUNTS,
            how="left",
            # TODO: Add in support for POS-sensitive filtering
            on=["CHROM", "POS", "ID", "REF", "ALT"],
            suffixes=("_CLASHING", ""),
        )
        _logger.info(
            "The %s populations allele count reports have been filed and processed.",
            POPULATION,
        )
        _logger.debug(DATA)

    _logger.info("Completed merging population-level allele count reports.")

    # %%
    ############ SAVE COUNT DATA TO FILE ############
    #################################################
    # for cluster in CLUSTERS:
    #     for gene in GENES:
    DATA.to_csv(snakemake.output[0], index=False)
    _logger.info("Results written to file.")
    # save_or_append_to_excel(
    #     CONFIG["output-dir"], DATA[cluster][gene], cluster, gene, "Count", "replace"
    # )
# %%

except Exception as E:
    _logger.error(E)
