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
import traceback

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


def chunks_of_n_size(l: iter, n: int):
    """
    This function splits a list into a generator of chunks of size n of the original list.
    """
    return (l[index : index + n] for index in range(0, len(l)))


def overlap(listA: list, listB: list):
    """
    This function determines if there is any overlap between two lists.
    """
    return list(set(listA).intersection(listB))


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
        names=MULTIINDEX + ["QUAL", "FILTER"],
        dtype={
            "CHROM": "Int8",
            "POS": "Int64",
            "ID": "string",
            "REF": "string",
            "ALT": "string",
            "QUAL": "string",
            "FILTER": "string",
        },
    )[MULTIINDEX]
    _logger.info("Variant index data has been imported.")
    _logger.debug(DATA)

    # _logger.info("Creating temporary directory for sequential I/O operations.")
    # DATA.to_csv(join(snakemake.output.tmp, "batch_0_p*.csv"), index=False)
    # _logger.info("Created temporary directory for sequential I/O operations.")
    # %%
    ############ MERGE RAW COUNTS ############
    ##########################################
    # [REPEAT] for each allele-count report produced per-population in this cluster

    _logger.info(
        "Beginning population count report files from available files found: [%s].",
        snakemake.input.allele_counts,
    )
    # for _batch_index, _report_batch in enumerate(
    #     chunks_of_n_size(snakemake.input.allele_counts, 5)
    # ):
    for report in snakemake.input.allele_counts:
        # Import data
        # DATA = dd.read_csv(
        #     join(report, f"batch_{_batch_index}_p*.csv"),
        #     dtype={
        #         "CHROM": "Int8",
        #         "POS": "Int64",
        #         "ID": "string",
        #         "REF": "string",
        #         "ALT": "string",
        #         "REF_CT": "Int64",
        #         "ALT_CTS": "Int64",
        #         "OBS_CT": "Int64",
        #     },
        # )

        # [FOR] each report in the batch:
        # for report in _report_batch:

        # Identify reported population.
        POPULATION = search("allele_count.([a-zA-Z]+).acount", report).group(1)
        _logger.info(
            "The %s population allele count report has been identified for analysis.",
            POPULATION,
        )

        # [IMPORT] reported populations allele count report
        REPORT = (
            read_csv(
                report,
                on_bad_lines="warn",
                sep="\t",  # This imports the columns under the following names. Length equals number of columns
                dtype={
                    "#CHROM": "Int8",
                    "POS": "Int64",
                    "ID": "string",
                    "REF": "string",
                    "ALT": "string",
                    "REF_CT": "Int64",
                    "ALT_CTS": "Int64",
                    "OBS_CT": "Int64",
                },
            )
            .rename(
                columns={
                    "#CHROM": "CHROM",
                    "ALT_CTS": f"{POPULATION}_ac",
                    "OBS_CT": f"{POPULATION}_tc",
                },
            )
            .drop(columns=["REF_CT"])
        )
        _logger.info(
            "The %s populations allele count report has been imported.",
            POPULATION,
        )
        _logger.info(REPORT)

        # [MERGE(left)] allele-count report into base `Dataframe`
        DATA = DATA.merge(
            REPORT,
            how="left",
            # TODO: Add in support for POS-sensitive filtering
            on=MULTIINDEX,
            suffixes=("_CLASHING", ""),
        )
        _logger.info(
            "The %s populations allele count reports have been imported.",
            POPULATION,
        )
        _logger.info(DATA)

    # [END(for)] a single report batch

    # Save this batch back to file and rinse and repeat until batches are all processed.
    # _logger.info("Writing results to file.")
    # DATA.to_csv(
    #     join(snakemake.output.tmp, f"batch_{_batch_index+1}_p*.csv"), index=False
    # )
    # _logger.info("Results written to file.")

    # _logger.debug("Writing results to file.")
    # DATA.to_csv(snakemake.output[0], index=False, single_file=True)
    # _logger.info("Results written to file.")

    _logger.info("Completed merging population-level allele count reports.")
    # %%
    ############ SAVE COUNT DATA TO FILE ############
    #################################################
    # for cluster in CLUSTERS:
    #     for gene in GENES:
    _logger.info("Writing results to file.")
    DATA.to_csv(snakemake.output["allele_counts"], index=False)
    _logger.info("Results written to file.")
    # save_or_append_to_excel(
    #     config["output"], DATA[cluster][gene], cluster, gene, "Count", "replace"
    # )
# %%

except Exception as E:
    _logger.error(E)
    _logger.error(traceback.format_exc())
    exit(1)
