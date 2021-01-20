# %%
import pandas as pd

#%%
# Declare Constants and Functions:
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']


def freq(alt: int, total: int) -> int:
    """Calculate percentage frequency for alleles based on count data.

    Args:
        alt (int): Number of observations of the alternate allele.
        total (int): Total number of observations of any allele.

    Returns:
        int: decimal percentage allele frequency.
    """
    if alt != 0 and total != 0:
        return row['ALT_CTS']/row['OBS_CT']
    if (alt == 0 and total != 0) or (alt == 0 and total == 0):
        return 0
    if alt != 0 and total == 0:
        return 999
#%%
# Load the Supplementary Table:
supplementary = dict()
for gene in genes:
    supplementary[gene] = pd.read_excel("../final/{}.xlsx".format(gene), sheet_name='VEP')[['ID', 'POS', 'REF', 'ALT']]
#%%

frequency_data = dict()
fishers_data = dict()

for gene in genes:
    fishers_data[gene] = supplementary[gene][['ID', 'POS', 'REF', 'ALT']]
    for pop in populations:
        supplementary[gene][pop] = 0
        fishers_data[gene]['{pop}_ac'.format(pop=pop)] = 0
        fishers_data[gene]['{pop}_tc'.format(pop=pop)] = 0
        frequency_data = pd.read_csv("../final/SUPER/ALL_{gene}.{pop}.acount".format(gene=gene, pop=pop), delimiter="\t").rename(columns={'#CHROM': 'CHROM'})
        for index, row in frequency_data.iterrows():
            supplementary[gene].loc[supplementary[gene]['ID'] == row['ID'], pop] = freq(row['ALT_CTS'], row['OBS_CT'])
            fishers_data[gene].loc[supplementary[gene]['ID'] == row['ID'], '{pop}_ac'.format(pop=pop)] = row['ALT_CTS']
            fishers_data[gene].loc[supplementary[gene]['ID'] == row['ID'], '{pop}_tc'.format(pop=pop)] = row['OBS_CT']

# %%
# Save the resulting dataframe back to its excel file:
for gene in genes:
    supplementary[gene].to_csv("../final/{}_Freq.csv".format(gene), index=False, sep="\t")

# %%
