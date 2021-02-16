# %%

# Import dependancies
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# %%

# Declare Constants and Functions
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
refPop = 'AFR'
compPop = ['AMR', 'EUR', 'EAS', 'SAS']

def Fishers (input: dict, refPop: str, compPop: list):
    """Runs a row-wise Fisher's Exact Test between the two listed populations


    Parameters:
    input (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data. 
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """
    for key, dataset in input.items():
        dataset['OR'] = dataset['Count'][['ID', 'POS', 'REF', 'ALT']]
        dataset['P'] = dataset['Count'][['ID', 'POS', 'REF', 'ALT']]
        for pop in compPop:
            for key2, direction in {"L": 'less', "G": 'greater', "T": "two-sided"}.items():
                oLabel = "{refPop}_{tail}_{pop}".format(refPop=refPop, tail=key2, pop=pop)
                pLabel = "{refPop}_{tail}_{pop}".format(refPop=refPop, tail=key2, pop=pop)
                dataset['OR'][oLabel] = None
                dataset['P'][pLabel] = None
                for index, row in dataset['Count'].iterrows():
                    contingency = [
                        [row["{}_ac".format(refPop)], row["{}_tc".format(refPop)]-row["{}_ac".format(refPop)]], 
                        [row["{}_ac".format(pop)], row["{}_tc".format(pop)]-row["{}_ac".format(pop)]]
                        ]
                    oVal, pVal = fisher_exact(contingency, alternative=direction)
                    dataset['P'].loc[index, pLabel] = pVal
                    dataset['OR'].loc[index, oLabel] = oVal
                dataset['P'][pLabel] = multipletests(dataset['P'][pLabel], method="bonferroni")[1]
        columnsToDrop = list()
        # for pop in populations:
        #     columnsToDrop.append("{}_ac".format(pop))
        #     columnsToDrop.append("{}_tc".format(pop))
        # dataset.drop(columns=columnsToDrop ,inplace=True)
# %%

# Load the Supplementary Table
supplementary = dict()
for gene in genes:
    supplementary[gene] = dict()
    supplementary[gene]['Count'] = pd.read_csv("../final/Supplementary Table/{}_Count.csv".format(gene), sep="\t")
    supplementary[gene]['OR'] = dict()
    supplementary[gene]['P'] = dict()

# %%

# Run Fisher's Exact test
Fishers(supplementary, refPop, compPop)

# %%

# Save Fishers Data to CSV
for gene in genes:
    supplementary[gene]["P"].to_csv("../final/Supplementary Table/{}_FishersP.csv".format(gene), sep='\t', index=False)
    supplementary[gene]['OR'].to_csv("../final/Supplementary Table/{}_FishersOR.csv".format(gene), sep='\t', index=False)
# %%
