#!/usr/bin/env python
"""
A Python script designed to run Frequency calculations and save the results.
"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from os.path import join

import pandas as pd
from pandas import read_csv, read_excel

from common.common import calculate_frequency
from common.fishers_exact_test import fishers_exact_test

############ IMPORT DATA ############
#####################################

# [IMPORT] sample-level data:
SAMPLE_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = pd.read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = pd.read_csv(join("..", "..", "input", "transcripts.csv"))


# %%
############ SET CONSTANT VARIABLES ############
################################################

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

# [SET] variant effect prediction metadata
ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# %%
# [SAVE] Frequency analysis to excel
# frequency_data = dict()
# fishers_data = dict()
supplementary = dict()

for cluster in CLUSTERS:
    # frequency_data[cluster] = dict()
    # fishers_data[cluster] = dict()
    supplementary[cluster] = dict()

    
    for gene in GENES:
        DATA_TEMPLATE = read_excel(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}.xlsx",
            ),
            sheet_name="Variant-EFfect-Prediction",
        )[["ID", "POS", "REF", "ALT"]]





        fishers_data[cluster][gene] = read_excel(
            join(
                "..",
                "..",
                "results",
                "FINAL",
                f"{cluster}-{gene}.xlsx",
            ),
            sheet_name="Variant-EFfect-Prediction",
        )
        for population_group in SAMPLE_DATA[cluster].unique():
            supplementary[cluster][gene][population_group] = 0
            supplementary[cluster][gene][f"{population_group}_ac"] = 0
            supplementary[cluster][gene][f"{population_group}_tc"] = 0
            frequency_data[cluster][gene] = read_csv(
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

            for index, row in frequency_data[cluster][gene].iterrows():
                supplementary[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"], population_group
                ] = calculate_frequency(row["ALT_CTS"], row["OBS_CT"])
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    f"{population_group}_ac",
                ] = row["ALT_CTS"]
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    f"{population_group}_tc",
                ] = row["OBS_CT"]


# Save the resulting dataframe back to its excel file:
for cluster in CLUSTERS:
    for gene in GENES:
        supplementary[cluster][gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_Freq.csv",
            ),
            index=False,
            sep="\t",
        )
        fishers_data[cluster][gene].to_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_Count.csv",
            ),
            index=False,
            sep="\t",
        )


# FISHERS EXACT


# Load the Supplementary Table
supplementary = dict()

for cluster in CLUSTERS:
    supplementary[cluster] = dict()
    for gene in GENES:
        supplementary[cluster][gene] = dict()
        supplementary[cluster][gene]["Count"] = pd.read_csv(
            join(
                "..",
                "..",
                "results",
                "Supplementary Table",
                cluster,
                f"{gene}_Count.csv",
            ),
            sep="\t",
        )
        supplementary[cluster][gene]["OR"] = dict()
        supplementary[cluster][gene]["P"] = dict()


# Run Fisher's Exact test
fishers_exact_test(supplementary["SUPER"], REFERENCE_POPULATION, COMPARISON_POPULATIONS)


# Save Fishers Data to CSV
for gene in GENES:
    supplementary["SUPER"][gene]["P"].to_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            "SUPER",
            f"{gene}_FishersP.csv",
        ),
        sep="\t",
        index=False,
    )
    supplementary["SUPER"][gene]["OR"].to_csv(
        join(
            "..",
            "..",
            "results",
            "Supplementary Table",
            "SUPER",
            f"{gene}_FishersOR.csv",
        ),
        sep="\t",
        index=False,
    )
