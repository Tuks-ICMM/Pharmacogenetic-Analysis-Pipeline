# %%
import pandas as pd
import seaborn as sns
# %%


# Declare Constants and Functions:
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
tests = ['VEP', 'Freq', 'Fishers']

# %%

data = dict()

for gene in genes:
    data[gene] = dict()
    for test in tests:
        data[gene][test] = pd.read_excel("../final/{}.xlsx".format(gene), sheet_name=test)

# %%

haplotypes = pd.read_excel("../final/Haplotypes.xlsx", sheet_name='table')
# %%

filteredData = dict()

for gene in genes:
    filteredData[gene] = dict()
    for test in tests:
        filteredData[gene][test] = data[gene][test].merge(haplotypes, left_on='ID', right_on='Core SNP')
# %%
