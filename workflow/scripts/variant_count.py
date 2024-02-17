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

from os.path import join

from common.common import read_vcf, save_or_append_to_excel
from pandas import merge, read_csv

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"


# %%
############ IMPORT DATA & SET CONSTANTS ############
#####################################################

# [IMPORT] sample-level data:
SAMPLE_DATA = read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = read_csv(join("..", "..", "input", "transcripts.csv"))

# [SET] cluster names:
CLUSTERS = set(SAMPLE_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# [SET] standard multi-index Columns
MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

# [IMPORT] variant count data:
DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        DATA[cluster][gene] = read_vcf(
            join("..", "..", "results", "FINAL", f"{cluster}", f"ALL_{gene}.vcf.gz")
        )[MULTIINDEX]
        # Rename the chromosome notation
        # DATA[cluster][gene]["CHROM"] = DATA[cluster][gene]["CHROM"].str.extract(
        #     "chr([1-9]{1,2}|[XY])"
        # )

# %%
############ MERGE RAW COUNTS ############
##########################################
for cluster in CLUSTERS:
    for gene in GENES:
        for population in SAMPLE_DATA[cluster].unique():
            TEMP_DATA = read_csv(
                join(
                    "..",
                    "..",
                    "results",
                    "FINAL",
                    cluster,
                    f"ALL_{gene}.{population}.acount",
                ),
                sep="\t",
                dtype={
                    "#CHROM": str,
                    "ID": str,
                    "REF": str,
                    "ALT": str,
                    "ALT_CTS": int,
                    "OBS_CT": int,
                },
            )[["#CHROM", "ID", "REF", "ALT", "ALT_CTS", "OBS_CT"]].rename(
                columns={
                    "#CHROM": "CHROM",
                    "ALT_CTS": f"{population}_ac",
                    "OBS_CT": f"{population}_tc",
                }
            )
            DATA[cluster][gene] = merge(
                DATA[cluster][gene],
                TEMP_DATA,
                how="inner",
                # ToDo: Add in support for POS-sensitive filtering
                on=["CHROM", "ID", "REF", "ALT"],
                suffixes=("_CLASHING", ""),
            )
# %%
############ SAVE COUNT DATA TO FILE ############
#################################################
for cluster in CLUSTERS:
    for gene in GENES:
        save_or_append_to_excel(DATA[cluster][gene], cluster, gene, "Count", "replace")
# %%
