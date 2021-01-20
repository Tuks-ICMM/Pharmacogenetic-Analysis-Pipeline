# %%
import pandas as pd
from scipy.stats import fisher_exact

# %%
# Declare Constants and Functions:
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
refPop = 'AFR'
compPop = ['AMR', 'EUR', 'EAS', 'SAS']

def Fishers (dict: dict, refPop: str, compPop: list):
    """Runs a row-wise Fisher's Exact Test between the two listed populations


    Parameters:
    dataframe (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data. 
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """
    for key, dataset in dict.items():
        for pop in compPop:
            oLabel = str(refPop) + "_OR_" + str(pop)
            pLabel = str(refPop) + "_P_" + str(pop)
            dataset[pLabel] = None
            dataset[oLabel] = None
            for index, row in dataset.iterrows():
                contingency = [
                    [row["{}_ac".format(refPop)], row["{}_tc".format(refPop)]-row["{}_ac".format(refPop)]], 
                    [row["{}_ac".format(pop)], row["{}_tc".format(pop)]-row["{}_ac".format(pop)]]
                    ]
                oVal, pVal = fisher_exact(contingency)
                dataset.loc[index, pLabel] = pVal
                dataset.loc[index, oLabel] = oVal
# %%
# Load the Supplementary Table:
supplementary = dict()
for gene in genes:
    supplementary[gene] = pd.read_excel("../final/{}.xlsx".format(gene), sheet_name='Fishers Exact')

# %%
# Run Fisher's Exact test on the datasets:
Fishers(supplementary, refPop, compPop)

# %%
for gene in genes:
    supplementary[gene].to_csv("../final/{}_Fishers.csv".format(gene), sep='\t', index=False)
# %%
