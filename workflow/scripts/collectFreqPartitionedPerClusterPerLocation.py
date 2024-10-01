#!/usr/bin/env python
"""
A Python script designed to run Frequency calculations and save the results to an Excel file.

[NEEDS] variant_count.py
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import logging

from common.common import calculate_frequency
from pandas import read_csv, cut

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

    SAMPLES = read_csv(snakemake.input.psam, sep="\t").rename(
        columns={"#IID": "sample_name"}
    )

    # [BUILD] list of populations to calculate frequencies for
    POPULATIONS = [
        column.replace("_tc", "") for column in ALLELE_COUNTS if column.endswith("_tc")
    ]
    _logger.info(
        "The following populations have been identified for analysis: %s", POPULATIONS
    )

    # [BUILD] Allele-frequency report from copy of allele-count report
    ALLELE_FREQUENCIES = ALLELE_COUNTS.copy(deep=True)

    # CREDIT: https://stackoverflow.com/questions/58454612/binning-a-column-in-a-dataframe-into-10-percentiles
    POPULATION_REPRESENTATION_LABELS = ["0-10%"] + [
        f"{bin_start+1}-{bin_start+10}%" for bin_start in range(10, 100, 10)
    ]

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
        _logger.info(
            "Calculating total number of samples provided for %s population", population
        )
        SAMPLE_NUMBER = SAMPLES.loc[
            SAMPLES[snakemake.wildcards.cluster] == population, "sample_name"
        ].count()
        _logger.info(
            "Calculated total number of samples provided for %s population", population
        )

        _logger.info(
            "Calculating population-specific representation of all loci observations given the sample annotations provided for %s population.",
            population,
        )
        ALLELE_FREQUENCIES[f"{population}_representation"] = cut(
            (ALLELE_FREQUENCIES[f"{population}_tc"] / SAMPLE_NUMBER).to_list(),
            10,
            labels=POPULATION_REPRESENTATION_LABELS,
        )
        _logger.info(
            "Calculated population-specific representation of all loci observations given the sample annotations provided for %s population.",
            population,
        )
    _logger.info(
        "Calculating list of unique populations for cluster %s.",
        snakemake.wildcards.cluster,
    )
    ALL_POPULATIONS = SAMPLES[snakemake.wildcards.cluster].unique()
    _logger.info(
        "Calculated list of unique populations for cluster %s. Unique populations identified: %s",
        snakemake.wildcards.cluster,
        ALL_POPULATIONS,
    )

    _logger.info("Calculating total allele counts across all populations observed.")
    # TODO: Import and incorporate plink AllFreq rule.
    ALLELE_FREQUENCIES["Total_ac"] = ALLELE_FREQUENCIES[
        [f"{population}_tc" for population in ALL_POPULATIONS]
    ].sum(
        axis=1
    )  # sum over column
    _logger.info("Calculated total allele counts across all populations observed.")

    _logger.info(
        "Calculating total representation of all loci given the sample annotations provided. The following populations were identified: %s ",
        ALL_POPULATIONS,
    )
    ALLELE_FREQUENCIES[f"Total_representation"] = cut(
        (ALLELE_FREQUENCIES["Total_ac"] / SAMPLES["sample_name"].count()),
        10,
        labels=POPULATION_REPRESENTATION_LABELS,
    )
    _logger.info(
        "Calculated total representation of all loci given the sample annotations provided."
    )

    for population in POPULATIONS:
        _logger.info(
            "The Allele-Count data for %s is being removed from the output.", population
        )
        # [DELETE] the allele-count columns for this population (memory-management)
        ALLELE_FREQUENCIES.drop(
            [f"{population}_ac", f"{population}_tc"], axis=1, inplace=True
        )
        _logger.info(
            "The Allele-Count data for %s has been removed from the output.", population
        )

    ############ SAVE TO OUTPUT ############
    ########################################
    _logger.info("Writing results to file.")
    ALLELE_FREQUENCIES.to_csv(snakemake.output["file"], index=False)
    _logger.info("Results written to file.")

except Exception as E:
    _logger.error(E)
