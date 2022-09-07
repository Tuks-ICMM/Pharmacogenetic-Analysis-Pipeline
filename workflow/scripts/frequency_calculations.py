#!/usr/bin/env python
"""
A Python script designed to run Frequency calculations and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import exists, join

import pandas as pd
from numpy import nan
from pandas import read_csv, read_excel

from common.common import calculate_frequency, read_vcf, save_or_append_to_excel

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
SAMPLE_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))


# [SET] cluster names:
CLUSTERS = set(SAMPLE_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")


# [IMPORT] gene-level data:
GENE_DATA = pd.read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = pd.read_csv(join("..", "..", "input", "transcripts.csv"))

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [IMPORT] Data and initial data-templates with variant IDs and population columns:
DATA = dict()
for cluster in CLUSTERS:
    DATA[cluster] = dict()
    for gene in GENES:
        # [1] Import the data in
        DATA[cluster][gene] = read_vcf(
            join("..", "..", "results", "FINAL", f"ALL_{gene}.vcf.gz")
        )[["CHROM", "POS", "ID", "REF", "ALT"]]

        # Extract the Chromosome code only:
        DATA[cluster][gene]["CHROM"] = DATA[cluster][gene]["CHROM"].str.extract(
            "chr([1-9]{1,2}|[XY])"
        )

        # Set empty columns for our incomming population groupings:
        DATA[cluster][gene][SAMPLE_DATA[cluster].unique()] = nan

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# [SET] variant effect prediction metadata
ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


# %%
############ CALCULATE VARIANT FREQUENCIES ############
#######################################################

for cluster in CLUSTERS:
    for gene in GENES:
        for population_group in SAMPLE_DATA[cluster].unique():
            DATA[cluster][gene][population_group] = 0

            # Move this to import section
            frequency_data = read_csv(
                join(
                    "..",
                    "..",
                    "results",
                    "FINAL",
                    cluster,
                    f"ALL_{gene}.{population_group}.acount",
                ),
                delimiter="\t",
            ).rename(columns={"#CHROM": "CHROM"})

            # Use the .loc function to inject the frequency data for each variant:
            for index, row in frequency_data.iterrows():
                DATA[cluster][gene].loc[
                    DATA[cluster][gene]["ID"] == row["ID"], population_group
                ] = calculate_frequency(row["ALT_CTS"], row["OBS_CT"])

# %%
############ SAVE TO EXCEL ############
#######################################
for cluster in CLUSTERS:
    for gene in GENES:
        save_or_append_to_excel(DATA[cluster][gene], cluster, gene, "Freq", "replace")
# %%
