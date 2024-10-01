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
from json import loads
import sys
from typing import Generator

import pandas as pd

# from common.condel_score import condel_weighted_score
from pandas import DataFrame, Series, read_csv

from scripts.common.variant_effect_prediction_parsers import (
    collectConsequenceTypes,
    collectFeatureTypes,
    extractTranscriptConsequenceValue,
)

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

    # [IMPORT] variant count data:
    _logger.info("Dataset identified for imported.")
    DATA = read_csv(
        snakemake.input.vep_results,
        on_bad_lines="warn",
    )
    _logger.info("Dataset has been successfully imported.")
    DATA.rename(columns={"#CHROM": "CHROM"})
    _logger.info(f"Dataset columns have been renamed. New columns are: {DATA.keys().tolist()}")
    DATA["result"] = DATA["result"].apply(lambda vep_response: loads(vep_response))
    _logger.info("Dataset results column has been parsed into JSON-format.")
    DATA.set_index(MULTIINDEX, inplace=True)
    _logger.info("Dataset has been re-indexed using global genomic multi-index.")

    # [EXTRACT] the VEP values.
    _logger.info("Extracting Feature_type.")
    DATA["Feature_type"] = DATA["result"].apply(lambda row: collectFeatureTypes(row))
    _logger.info("Extracted Feature_type.")

    _logger.info("Extracting Consequence_type.")
    DATA["Consequence_type"] = DATA["result"].apply(
        lambda row: collectConsequenceTypes(row)
    )
    _logger.info("Extracted Consequence_type.")

    _logger.info("Extracting Biotype.")
    DATA["Biotype"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "biotype", __name__)
    )
    _logger.info("Extracted Biotype.")

    _logger.info("Extracting CADD.")
    DATA["CADD"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "cadd_raw", __name__)
    )
    _logger.info("Extracted CADD.")

    _logger.info("Extracting CADD_PHRED.")
    DATA["CADD_PHRED"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "cadd_phred", __name__)
    )
    _logger.info("Extracted CADD_PHRED.")

    _logger.info("Extracting Gene_symbol.")
    DATA["Gene_symbol"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(
            result, "gene_symbol", __name__
        )
    )
    _logger.info("Extracted Gene_symbol.")

    _logger.info("Extracting Gene_symbol_source.")
    DATA["Gene_symbol_source"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(
            result, "gene_symbol_source", __name__
        )
    )
    _logger.info("Extracted Gene_symbol_source.")

    _logger.info("Extracting Impact.")
    DATA["Impact"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "impact", __name__)
    )
    _logger.info("Extracted Impact.")

    _logger.info("Extracting HGVS.")
    DATA["HGVSC"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "hgvsc", __name__)
    )
    _logger.info("Extracted HGVS.")

    _logger.info("Extracting Gene_ID.")
    DATA["Gene_ID"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "gene_id", __name__)
    )
    _logger.info("Extracted Gene_ID.")

    _logger.info("Extracting Transcript_ID.")
    DATA["Transcript_ID"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(
            result, "transcript_id", __name__
        )
    )
    _logger.info("Extracted Transcript_ID.")

    _logger.info("Extracting Exon.")
    DATA["Exon"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "exon", __name__)
    )
    _logger.info("Extracted Exon.")

    _logger.info("Extracting isCanonical.")
    DATA["isCanonical"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "canonical", __name__)
    )
    _logger.info("Extracted isCanonical.")

    _logger.info("Extracting Variant_allele.")
    DATA["Variant_allele"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(
            result, "variant_allele", __name__
        )
    )
    _logger.info("Extracted Variant_allele.")

    _logger.info("Dropping unneeded results column.")
    DATA.drop(columns=["result"], inplace=True)
    _logger.info("Extracted Variant_allele.")

    _logger.info(
        "Saving results to output file: '%s'.", snakemake.output.cleaned_vep_results
    )
    DATA.to_csv(snakemake.output.cleaned_vep_results)
    _logger.info(
        "Saved results to output file: '%s'.", snakemake.output.cleaned_vep_results
    )


# %%

except Exception as E:
    _logger.error(E)
    sys.exit()
