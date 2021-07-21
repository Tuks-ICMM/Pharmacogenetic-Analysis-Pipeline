# %%

# Import dependancies
import pandas as pd
import json
from os.path import join

#%%

# Declare Constants and Functions:
with open(join("..", "..", "config", "config.json")) as f:
    config = json.load(f)
genes = config["locations"]
clusters = config["cluster"]["clusters"]
popsFile = pd.read_excel(join("..", "..", "config", "Clusters.xlsx"))
populations = ["AFR", "AMR", "EUR", "EAS", "SAS"]


def freq(alt: int, total: int) -> int:
    """Calculate percentage frequency for alleles based on count data.

    Args:
        alt (int): Number of observations of the alternate allele.
        total (int): Total number of observations of any allele.

    Returns:
        int: decimal percentage allele frequency.
    """
    if alt != 0 and total != 0:
        return row["ALT_CTS"] / row["OBS_CT"]
    if (alt == 0 and total != 0) or (alt == 0 and total == 0):
        return 0
    if alt != 0 and total == 0:
        return 999


#%%

# Load the Supplementary Table:
supplementary = dict()

for cluster in clusters:
    supplementary[cluster] = dict()
    for gene in genes:
        supplementary[cluster][gene] = pd.read_csv(
            "../final/Supplementary Table/{cluster}/{gene}_VEP.csv".format(
                cluster=cluster, gene=gene
            ),
            sep="\t",
        )[["ID", "POS", "REF", "ALT"]]

#%%

# Calculate frequencies
frequency_data = dict()
fishers_data = dict()


for cluster in clusters:
    frequency_data[cluster] = dict()
    fishers_data[cluster] = dict()
    for gene in genes:
        fishers_data[cluster][gene] = supplementary[cluster][gene][
            ["ID", "POS", "REF", "ALT"]
        ]
        for pop in popsFile[cluster].unique():
            # supplementary[cluster][gene][pop] = 0
            fishers_data[cluster][gene]["{pop}_ac".format(pop=pop)] = 0
            fishers_data[cluster][gene]["{pop}_tc".format(pop=pop)] = 0
            frequency_data[cluster][gene] = pd.read_csv(
                "../final/{cluster}/ALL_{gene}.{pop}.acount".format(
                    cluster=cluster, gene=gene, pop=pop
                ),
                delimiter="\t",
            ).rename(columns={"#CHROM": "CHROM"})
            for index, row in frequency_data[cluster][gene].iterrows():
                supplementary[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"], pop
                ] = freq(row["ALT_CTS"], row["OBS_CT"])
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    "{pop}_ac".format(pop=pop),
                ] = row["ALT_CTS"]
                fishers_data[cluster][gene].loc[
                    supplementary[cluster][gene]["ID"] == row["ID"],
                    "{pop}_tc".format(pop=pop),
                ] = row["OBS_CT"]

# %%

# Save the resulting dataframe back to its excel file:
for cluster in clusters:
    for gene in genes:
        supplementary[cluster][gene].to_csv(
            "../final/Supplementary Table/{cluster}/{gene}_Freq.csv".format(
                cluster=cluster, gene=gene
            ),
            index=False,
            sep="\t",
        )
        fishers_data[cluster][gene].to_csv(
            "../final/Supplementary Table/{cluster}/{gene}_Count.csv".format(
                cluster=cluster, gene=gene
            ),
            index=False,
            sep="\t",
        )


# %%
