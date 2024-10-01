#!/usr/bin/env python
"""
A Python script designed to run Variant Effect Prediction calculations
and API calls as well as save the results to file.

[NEEDS] WORKFLOW_RUNTIME

"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

import logging
from json import dumps
from typing import Generator
from time import sleep

import pandas as pd
from common.common import chunk, generate_notation, generate_params

# from common.condel_score import condel_weighted_score
from pandas import DataFrame, Series, read_csv
from requests import post

# from workflow.scripts.entities.VariantConsequence import VariantConsequenceResult

# [SET] Pandas Chained-Assignment off to prevent writing changes to transient DF copies:
pd.set_option("chained_assignment", None)

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Fatima Barmania",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

LOGGER = logging.getLogger("variant_count.py")
LOGGER.setLevel(logging.DEBUG)

LOG_FILE = logging.FileHandler(snakemake.log[0])
LOG_FILE.setLevel(logging.DEBUG)
LOGGER.addHandler(LOG_FILE)


def generate_notation(index: Series, gene_strand: str) -> str:
    """
    Parses a row and returns the E! Ensemble query notation for a provided mutation.

    Args:
        row (Series): A row (generated from .iterrows() method) for which notation is needed.
    Returns:
        str: HGVS notation for the given variant row.
    """
    # TODO: Review this algorithm.
    CHROM, POS, ID, REF, ALT = index
    if (len(REF) > len(ALT)) | (len(REF) == len(ALT)):
        stop_coordinates = int(POS) + (len(REF) - 1)
        return f"{CHROM}:{POS}-{stop_coordinates}:1/{ALT}"  # TODO: Check strand formatting here. Appears to be incorectly hard-coded
    elif len(REF) < len(ALT):
        stop_coordinates = int(POS) + len(REF)
        return f"{CHROM}:{POS}-{stop_coordinates}:{gene_strand}/{ALT}"


def chunk(input_data: DataFrame, size: int) -> Generator:
    """
    Render a generator to yield chunks of data from the original file.

    Args:
        input_data (DataFrame): The dataset to divide into chunks.
        size (int): The maximum size of the resulting chunks.
    Returns:
        [type]: [description]
    Yields:
        Generator: [description]
    """
    return (
        input_data[pos : pos + size]["query"] for pos in range(0, len(input_data), size)
    )


# %%
try:
    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [SET] variant effect prediction metadata
    ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
    HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

    # [IMPORT] variant count data:
    DATA = read_csv(
        snakemake.input.pvar,
        on_bad_lines="warn",
        comment="#",
        sep="\t",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
    )[MULTIINDEX]

    LOGGER.debug("Dataset has been imported")

    DATA["query"] = ""
    DATA.set_index(MULTIINDEX, inplace=True)
    LOGGER.debug(f"Multiindex has been set: {DATA.index}")

    ############ GENERATE HGVS NOTATION ############
    ################################################

    DATA["query"] = DATA.apply(
        lambda row: generate_notation(row.name, snakemake.params.strand), axis="columns"
    )
    LOGGER.debug("E! Ensemble query notation has been compiled")

    # [FOR] each batch of 200 variants,
    for payload in chunk(DATA, 200):
        LOGGER.debug(
            f"Attempting to query payload: {dumps(dict(variants=payload.tolist()))}"
        )
        # [QUERY] the batch and store the E! Ensemble API response
        response = post(
            ENDPOINT,
            headers=HEADERS,
            params=generate_params(),
            data=dumps(dict(variants=payload.tolist())),
        )
        LOGGER.debug("Response received")

        # [IF] response received was valid:
        if response.ok:
            # [FOR] record
            for record in response.json():
                DATA.loc[DATA["query"] == record["input"], "result"] = dumps(record)
            LOGGER.debug("API batch has completed successfully")
        else:
            LOGGER.warn(f"API call did not come back as ok: {response}")
        sleep(2)
    DATA.to_csv(snakemake.output[0])

# %%

except Exception as E:
    LOGGER.error(E)
