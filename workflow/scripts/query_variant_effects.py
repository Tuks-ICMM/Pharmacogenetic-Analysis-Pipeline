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
import traceback
from typing import Generator
from time import sleep

import pandas as pd
from common.common import chunk, generate_notation

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

VEP_QUERY_PARAMETERS = {
        # https://www.ensembl.org/info/genome/compara/epo_anchors_info.html
        # Retrive the ancestral allele at this locus as per the EPO pipeline
        # "AncestralAllele": True,
        # Amino Acid conservation score (Blosum62 method)
        "Blosum62": True,
        # Request a CADD score for this variant
        "CADD": True,
        # Retrive a conservation score from the E! Ensembl Compara database
        # "Conservation": True,
        # https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/master/DisGeNET.pm
        # Retrieves a list of variant-disease PMID associations for meta-analysis
        # "DisGeNET": True,
        "Enformer": True,
        # Retrieve variant classification using evolutionary sequences
        # "EVE": True,
        # Retrive Gene-Ontology terms associated with the variant
        # "GO": True,
        # Retrive splice sites associated with this variant
        # "GeneSplicer": True,
        # Retrive phenotypic profiles for a variant sequence defined using human phenotype ontology terms
        # "Geno2MP": True,
        # https://www.ebi.ac.uk/intact/home
        # Retrive a list of molecular interactions involving this variant asper the IntAct database
        # "IntAct": True,
        # Retrive an indicator for Loss-of-function for the given variant.
        "LoF": True,
        # https://www.genomenon.com/mastermind/
        # Retrive a list of associated literature which cites the variant using the MasterMind database
        "Mastermind": True,
        # Retrive a score from the MaveDB database based on multiplex assay datasets
        # "MaveDB": True,
        # Retrive splice-site consensus predictions based on maximum entropy
        "MaxEntScan": True,
        # Predict if a variant allows nonsense-mediated mRNA decay
        # "NMD": True,
        # Retrives phenotype records that overlap
        # "Phenotypes": True,
        # Retrives pre-calculated SpliceAI to predict splice junctions. I have selected 2 to pull MANE annotations.
        # "SpliceAI": 2,
        # Predicts impact of 5' UTR variants (New ORFs, etc)
        "UTRAnnotator": True,
        # Retrive APRIS isoform information for the given variant.
        # "appris": True,
        # Request that canonical transcripts be flagged
        "canonical": True,
        # Retrive a list of CCDS identifiers for recognized protein-coding regions
        # "ccds": True,
        # Retrive pathogenicity predictions for the variant from dbNSFP
        "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MetaSVM_score,MetaSVM_pred,Aloft_Fraction_transcripts_affected,transcript_match=1",
        # Retrive rpedictions for splice variants
        # "dbscSNV": True,
        # Request a list of overlapping protein domain names
        "domains": True,
        # Outputs only the most severe consequence per gene, using the criteria set by 'pick_order'
        "per_gene": True,
        # Selects single consequence record based on variant allele and gene combination
        "pick_allele_gene": True,
        # Pick one line or block of consequence data per variant, including transcript-specific columns.
        "pick": True,
        # Select the criteria order to use when selecting a single consequence using the 'pick' flag
        "pick_order": "canonical,mane_plus_clinical,mane_select,appris,tsl,biotype,ccds,rank,length",
        # Requests GA4GH Variation Representation Specification annotations
        # "ga4gh_vrs": True,
        # Request HGVS nomenclature
        "hgvs": True,
        # Request MANE Select annotations
        # "mane": True,
        # Retrive miRNA secondary structure annotations for this variant
        # "mirna": True,
        # Retrive predictions for destabilization effect of variant
        "mutfunc": True,
        # Retrive number fo affected intron and exon regions in transcript
        "numbers": True,
        # Retrive E! Ensembl protein identifiers
        # "protein": True,
        # retrive REVEL scores
        "REVEL": True,
        # Shift all variants that overlap ith transcripts as far possible in a 3' direction before predicting consequences
        # "shift_3prime": True,
        # Retrive transcript version numbers as well
        # "transcript_version": True,
        # Retrive transcript support level annotations regarding how well mRNA aligns across splice sites
        # "tsl": True,
        # Retrive accessions for the gene from three protein product databases
        # "uniprot": True,
        # Retrive variant class annotations based on sequence ontology
        "variant_class": True,
    }

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
        LOGGER.debug(f"Using the following parameters: {VEP_QUERY_PARAMETERS}")
        # [QUERY] the batch and store the E! Ensemble API response
        response = post(
            ENDPOINT,
            headers=HEADERS,
            params=VEP_QUERY_PARAMETERS,
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
    LOGGER.error(traceback.format_exc())
    exit(1)
