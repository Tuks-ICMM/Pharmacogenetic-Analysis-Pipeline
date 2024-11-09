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
import traceback
from typing import Generator

import pandas as pd

# from common.condel_score import condel_weighted_score
from pandas import DataFrame, Series, read_csv

from scripts.common.variant_effect_prediction_parsers import (
    collectConsequenceTypes,
    collectFeatureTypes,
    extractTranscriptConsequenceValue,
    extractMutFuncValues,
    extractUTRAnnotatorValues
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

    _logger.info("Extracting REVEL score.")
    DATA["REVEL_score"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "revel", __name__)
    )
    _logger.info("Extracted REVEL score.")

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

    _logger.info("Extracting SIFT 4G score.")
    DATA["SIFT_4G_score"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "sift4g_score", __name__)
    )
    _logger.info("Extracted SIFT 4G score.")


    _logger.info("Extracting SIFT 4G prediction.")
    DATA["SIFT_4G_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "sift4g_pred", __name__)
    )
    _logger.info("Extracted SIFT 4G prediction")


    _logger.info("Extracting SIFT 4G score.")
    DATA["SIFT_score"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "sift_score", __name__)
    )
    _logger.info("Extracted SIFT 4G score.")


    _logger.info("Extracting SIFT prediction.")
    DATA["SIFT_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "sift_pred", __name__)
    )
    _logger.info("Extracted SIFT prediction.")


    _logger.info("Extracting PolyPhen 4G score.")
    DATA["PolyPhen_score"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "polyphen_score", __name__)
    )
    _logger.info("Extracted PolyPhen 4G score.")


    _logger.info("Extracting PolyPhen prediction.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "polyphen_pred", __name__)
    )
    _logger.info("Extracted PolyPhen prediction.")

    _logger.info("Extracting FATHMM score.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "fathmm_score", __name__)
    )
    _logger.info("Extracted FATHMM score.")

    _logger.info("Extracting FATHMM prediction.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "fathmm_pred", __name__)
    )
    _logger.info("Extracted FATHMM prediction.")

    _logger.info("Extracting PROVEAN score.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "provean_score", __name__)
    )
    _logger.info("Extracted PROVEAN score.")

    _logger.info("Extracting PROVEAN prediction.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "provean_pred", __name__)
    )
    _logger.info("Extracted PROVEAN prediction.")

    _logger.info("Extracting MetaSVM score.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "metasvm_score", __name__)
    )
    _logger.info("Extracted MetaSVM score.")

    _logger.info("Extracting MetaSVM prediction.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "metasvm_pred", __name__)
    )
    _logger.info("Extracted MetaSVM prediction.")

    _logger.info("Extracting Aloft fraction of transcripts affected.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_fraction_transcripts_affected", __name__)
    )
    _logger.info("Extracted Aloft fraction of transcripts affected.")


    _logger.info("Extracting Aloft fraction of transcripts affected.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_fraction_transcripts_affected", __name__)
    )
    _logger.info("Extracted Aloft fraction of transcripts affected.")


    _logger.info("Extracting Aloft prob Tolerant.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_prob_tolerant", __name__)
    )
    _logger.info("Extracted Aloft prob Tolerant.")

    _logger.info("Extracting Aloft prob Recessive.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_prob_recessive", __name__)
    )
    _logger.info("Extracted Aloft prob Recessive.")

    _logger.info("Extracting Aloft prob Dominant.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_prob_dominant", __name__)
    )
    _logger.info("Extracted Aloft prob Dominant.")

    _logger.info("Extracting Aloft pred.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_pred", __name__)
    )
    _logger.info("Extracted Aloft pred.")

    _logger.info("Extracting Aloft Confidance.")
    DATA["PolyPhen_prediction"] = DATA["result"].apply(
        lambda result: extractTranscriptConsequenceValue(result, "aloft_confidance", __name__)
    )
    _logger.info("Extracted Aloft Confidance.")

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

    _logger.info("Extracting MutFunc dG_mt results.")
    DATA["mutfunc_mod_dG_mt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "dG_mt"
        )
    )
    _logger.info("Extracted MutFunc dG_mt results.")

    _logger.info("Extracting MutFunc ddG_sd results.")
    DATA["mutfunc_mod_ddG_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "ddG_sd"
        )
    )
    _logger.info("Extracted MutFunc ddG_sd results.")

    _logger.info("Extracting MutFunc dG_wt results.")
    DATA["mutfunc_mod_dG_wt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "dG_wt"
        )
    )
    _logger.info("Extracted MutFunc dG_wt results.")

    _logger.info("Extracting MutFunc dG_mt_sd results.")
    DATA["mutfunc_mod_dG_mt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "dG_mt_sd"
        )
    )
    _logger.info("Extracted MutFunc dG_mt_sd results.")

    _logger.info("Extracting MutFunc dG_wt_sd results.")
    DATA["mutfunc_mod_dG_wt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "dG_wt_sd"
        )
    )
    _logger.info("Extracted MutFunc dG_wt_sd results.")

    _logger.info("Extracting MutFunc ddG results.")
    DATA["mutfunc_mod_ddG"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "mod", "ddG"
        )
    )
    _logger.info("Extracted MutFunc ddG results.")

    _logger.info("Extracting MutFunc experimental dG_mt results.")
    DATA["mutfunc_exp_dG_mt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "dG_mt"
        )
    )
    _logger.info("Extracted MutFunc experimental dG_mt results.")

    _logger.info("Extracting MutFunc experimental ddG_sd results.")
    DATA["mutfunc_exp_ddG_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "ddG_sd"
        )
    )
    _logger.info("Extracted MutFunc experimental ddG_sd results.")

    _logger.info("Extracting MutFunc experimental dG_wt results.")
    DATA["mutfunc_exp_dG_wt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "dG_wt"
        )
    )
    _logger.info("Extracted MutFunc experimental dG_wt results.")

    _logger.info("Extracting MutFunc experimental dG_mt_sd results.")
    DATA["mutfunc_exp_dG_mt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "dG_mt_sd"
        )
    )
    _logger.info("Extracted MutFunc experimental dG_mt_sd results.")

    _logger.info("Extracting MutFunc experimental dG_wt_sd results.")
    DATA["mutfunc_exp_dG_wt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "dG_wt_sd"
        )
    )
    _logger.info("Extracted MutFunc experimental dG_wt_sd results.")

    _logger.info("Extracting MutFunc experimental ddG results.")
    DATA["mutfunc_exp_ddG"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "exp", "ddG"
        )
    )
    _logger.info("Extracted MutFunc experimental ddG results.")

    _logger.info("Extracting MutFunc intersection dG_mt results.")
    DATA["mutfunc_int_dG_mt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "dG_mt"
        )
    )
    _logger.info("Extracted MutFunc intersection dG_mt results.")

    _logger.info("Extracting MutFunc intersection ddG_sd results.")
    DATA["mutfunc_int_ddG_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "ddG_sd"
        )
    )
    _logger.info("Extracted MutFunc intersection ddG_sd results.")

    _logger.info("Extracting MutFunc intersection dG_wt results.")
    DATA["mutfunc_int_dG_wt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "dG_wt"
        )
    )
    _logger.info("Extracted MutFunc intersection dG_wt results.")

    _logger.info("Extracting MutFunc intersection dG_mt_sd results.")
    DATA["mutfunc_int_dG_mt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "dG_mt_sd"
        )
    )
    _logger.info("Extracted MutFunc intersection dG_mt_sd results.")

    _logger.info("Extracting MutFunc intersection dG_wt_sd results.")
    DATA["mutfunc_int_dG_wt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "dG_wt_sd"
        )
    )
    _logger.info("Extracted MutFunc intersection ddG dG_wt_sd results.")

    _logger.info("Extracting MutFunc intersection results.")
    DATA["mutfunc_int_ddG"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "int", "ddG"
        )
    )
    _logger.info("Extracted MutFunc intersection ddG results.")

    _logger.info("Extracting MutFunc motif dG_mt results.")
    DATA["mutfunc_motif_dG_mt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "dG_mt"
        )
    )
    _logger.info("Extracted MutFunc motif dG_mt results.")

    _logger.info("Extracting MutFunc motif ddG_sd results.")
    DATA["mutfunc_motif_ddG_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "ddG_sd"
        )
    )
    _logger.info("Extracted MutFunc motif ddG_sd results.")

    _logger.info("Extracting MutFunc motif dG_wt results.")
    DATA["mutfunc_motif_dG_wt"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "dG_wt"
        )
    )
    _logger.info("Extracted MutFunc motif dG_wt results.")

    _logger.info("Extracting MutFunc motif dG_mt_sd results.")
    DATA["mutfunc_motif_dG_mt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "dG_mt_sd"
        )
    )
    _logger.info("Extracted MutFunc motif dG_mt_sd results.")

    _logger.info("Extracting MutFunc motif dG_wt_sd results.")
    DATA["mutfunc_motif_dG_wt_sd"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "dG_wt_sd"
        )
    )
    _logger.info("Extracted MutFunc motif dG_wt_sd results.")

    _logger.info("Extracting MutFunc motif ddG results.")
    DATA["mutfunc_motif_ddG"] = DATA["result"].apply(
        lambda result: extractMutFuncValues(
            result, "motif", "ddG"
        )
    )
    _logger.info("Extracted MutFunc motif ddG results.")

    # START: UTRAnnotator
    _logger.info("Extracting UTRAnnotator 5UTR_consequence results.")
    DATA["UTRAnnotator_5UTR_consequence"] = DATA["result"].apply(
        lambda result: extractUTRAnnotatorValues(
            result, "5utr_consequence"
        )
    )
    _logger.info("Extracted UTRAnnotator 5UTR_consequence results.")

    _logger.info("Extracting UTRAnnotator 5UTR_annotation results.")
    DATA["UTRAnnotator_5UTR_annotation"] = DATA["result"].apply(
        lambda result: extractUTRAnnotatorValues(
            result, "5utr_annotation"
        )
    )
    _logger.info("Extracted UTRAnnotator 5UTR_annotation results.")

    _logger.info("Extracting UTRAnnotator Existing_uORFs results.")
    DATA["UTRAnnotator_Existing_uORFs"] = DATA["result"].apply(
        lambda result: extractUTRAnnotatorValues(
            result, "existing_uorfs"
        )
    )
    _logger.info("Extracted UTRAnnotator Existing_uORFs results.")

    _logger.info("Extracting UTRAnnotator 5UTR_annotation results.")
    DATA["UTRAnnotator_Existing_OutOfFrame_oORFs"] = DATA["result"].apply(
        lambda result: extractUTRAnnotatorValues(
            result, "existing_outofframe_oorfs"
        )
    )
    _logger.info("Extracted UTRAnnotator Existing_OutOfFrame_oORFs results.")

    _logger.info("Extracting UTRAnnotator Existing_InFrame_oORFs results.")
    DATA["UTRAnnotator_Existing_InFrame_oORFs"] = DATA["result"].apply(
        lambda result: extractUTRAnnotatorValues(
            result, "existing_inframe_oorfs"
        )
    )
    _logger.info("Extracted UTRAnnotator Existing_InFrame_oORFs results.")
    # END: UTRAnnotator

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
    _logger.error(traceback.format_exc())
    exit(1)
