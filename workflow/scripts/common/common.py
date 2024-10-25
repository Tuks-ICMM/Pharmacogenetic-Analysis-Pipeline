"""
A suite of general-purpose functions used in this pipeline.
"""

from gzip import open as gzip_open
from io import StringIO
from os import makedirs
from os.path import exists, join
from typing import Generator

from numpy import nan
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo
from pandas import DataFrame, ExcelWriter, Series, read_csv

# %%

#########################################################
############ DEFINE GENRAL-PURPOSE VARIABLES ############
#########################################################
MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

############################################################
############ DEFINE GENRAL-PURPOSE CALCULATIONS ############
############################################################
def merge(input_row: set) -> str:
    """Generate row value to save in df
    Args:
        value (set): The set you wish to convert
    Returns:
        str: The converted value, ready to store.
    """
    row_value_as_list = list(input_row)
    if len(row_value_as_list) == 1:
        return str(row_value_as_list[0])
    else:
        return " | ".join(input_row)


def chunk(input_data: DataFrame, size: int) -> Generator:
    """Renders a generator to yield chunks of data from the original file.
    Args:
        input_data (DataFrame): The dataset to divide into chunks.
        size (int): The maximum size of the resulting chunks.
    Returns:
        [type]: [description]
    Yields:
        Generator: [description]
    """
    return (input_data[pos : pos + size] for pos in range(0, len(input_data), size))


def generate_notation(input_row: Series, gene_strand: str) -> str:
    """
    Parses a row and returns the HGVS notation to query E! Ensembl.

    Args:
        row (Series): A row (generated from .iterrows() method) for which notation is needed.
    Returns:
        str: HGVS notation for the given variant row.
    """
    CHROM, POS, REF, ALT = input_row.name
    if (len(REF) > len(ALT)) | (len(REF) == len(ALT)):
        stop_coordinates = int(POS) + (len(REF) - 1)
        return f"{CHROM}:{POS}-{stop_coordinates}:1/{ALT}"
    elif len(input_row["REF"]) < len(input_row["ALT"]):
        stop_coordinates = int(input_row["POS"]) + len(input_row["REF"])
        return f"{CHROM}:{POS}-{stop_coordinates}:{gene_strand}/{ALT}"


def read_vcf(path: str) -> DataFrame:
    """
    Function to import a VCF file as a Pandas dataframe. Does not include Genotype columns.
    CREDIT: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

    Args:
        path (str): Path to the VCF file in question
    Returns:
        Dataframe: A Pandas Dataframe of the VCF files content.
    """
    with gzip_open(path, "rt") as file:
        lines = [l for l in file if not l.startswith("##")]
    return read_csv(
        StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})


def calculate_frequency(alternate_allele_count: int, total_count: int) -> int:
    """Calculate percentage frequency for alleles based on count data.
    Args:
        alt (int): Number of observations of the alternate allele.
        total (int): Total number of observations of any allele.
    Returns:
        int: decimal percentage allele frequency.
    """
    if alternate_allele_count != 0 and total_count == 0:
        raise Exception("The total_count cannot be zero.")
    if alternate_allele_count == 0 and total_count == 0:
        return 0
    return alternate_allele_count / total_count


def generate_params(transcription_ids: list = None, canonical=False) -> dict:
    """
    Returns generated parameters suitable for the E! Ensembl VEP tool API.

    Args:
        transcription_ids (list): A list of transcript IDs to filter for.
        canonical (bool, optional): Only consider canonical transcripts. Defaults to False.

    Returns:
        dict: A dictionary of parameters suitable for the VEP API.
    """
    params = {
        # https://www.ensembl.org/info/genome/compara/epo_anchors_info.html
        # Retrive the ancestral allele at this locus as per the EPO pipeline
        # "AncestralAllele": True,
        # Amino Acid conservation score (Blosum62 method)
        "Blosum62": True,
        # Request a CADD score for this variant
        "CADD": True,
        # Retrive a conservation score from teh E! Ensembl Compara database
        # "Conservation": True,
        # https://raw.githubusercontent.com/ensembl-variation/VEP_plugins/master/DisGeNET.pm
        # Retrieves a list of variant-disease PMID associations for meta-analysis
        # "DisGeNET": True,
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
        # "MaxEntScan": True,
        # Predict if a variant allows nonsense-mediated mRNA decay
        # "NMD": True,
        # Retrives phenotype records that overlap
        # "Phenotypes": True,
        # Retrives pre-calculated SpliceAI to predict splice junctions. I have selected 2 to pull MANE annotations.
        # "SpliceAI": 2,
        # Predicts impact of 5' UTR variants (New ORFs, etc)
        # "UTRAnnotator": True,
        # Retrive APRIS isoform information for the given variant.
        # "appris": True,
        # Request that canonical transcripts be flagged
        "canonical": canonical,
        # Retrive a list of CCDS identifiers for recognized protein-coding regions
        # "ccds": True,
        # Retrive pathogenicity predictions for the variant from dbNSFP
        "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MetaSVM_score,MetaSVM_pred,Aloft_Fraction_transcripts_affected,transcript_match=1",
        # Retrive rpedictions for splice variants
        # "dbscSNV": True,
        # Request a list of overlapping protein domain names
        # "domains": True,
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
        # "mutfunc": True,
        # Retrive number fo affected intron and exon regions in transcript
        "numbers": True,
        # Retrive E! Ensembl protein identifiers
        # "protein": True,
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
    if transcription_ids:
        params["transcription_ids"] = transcription_ids
    return params


def directory_exists(path: str):
    """Test weather or not a directory exists. If not, create it.
    Args:
        path (str): file path of the directory to test.
    """
    if not exists(path):
        makedirs(path)


# ToDo: Find better way to type hint exists_behaviour
def save_or_append_to_excel(
    path: list,
    data_to_save: DataFrame,
    cluster: str,
    gene: str,
    sheet_name: str,
    exists_behavior: str,
):
    """
    A small function that writes or appends to an Excel file depending on if it exists or not.

    Args:
        data_to_save (DataFrame): A Pandas DataFrame to save to Excel format.
        cluster (str): The cluster level being saved.
        gene (str): The gene being saved.
    """
    path = join(
        *path,
        "Excel",
        f"{cluster}-{gene}.xlsx",
    )
    if exists(path):
        mode = "a"
    else:
        mode = "w"
    kwargs = {"engine": "openpyxl", "mode": mode}
    if exists_behavior and mode == "a":
        kwargs["if_sheet_exists"] = exists_behavior
    with ExcelWriter(  # pylint: disable=abstract-class-instantiated
        path, **kwargs
    ) as writer:
        data_to_save.to_excel(
            writer,
            index=False,
            sheet_name=sheet_name,
        )
        worksheet = writer.sheets[sheet_name]
        (row_count, column_count) = data_to_save.shape

        # Define the table using the OpenPyExel engine:
        table = Table(
            ref=f"A1:{get_column_letter(column_count)}{row_count+1}",
            displayName=sheet_name,
        )

        # Add a nice style, because I like nice-looking tables!
        # Add a default style with striped rows and banded columns
        style = TableStyleInfo(
            name="TableStyleMedium9",
            showFirstColumn=False,
            showLastColumn=False,
            showRowStripes=True,
            showColumnStripes=True,
        )
        table.tableStyleInfo = style

        # Add an Excel table to the data:
        worksheet.add_table(table)


# %%
