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
        data[gene][test] = pd.read_excel("../final/CYP2A6.xlsx", sheet_name=test)
# %%

for gene in genes:
    data[gene]['Fishers']['AFR_rc'] = data[gene]['Fishers']['AFR_tc'] - data[gene]['Fishers']['AFR_ac']
    sns.displot(data[gene]['Fishers'], x='AFR_ac', y='AFR_rc', kind="kde", binwidth=(100, ))
# %%
