#%% 

%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import myvariant as mv
import numpy as np
import pyupset as pup

#%%

# Define datasets:

datasets = ['CYP2A6','CYP2B6', 'UGT2B7']
data = dict()

#%%


# for dataset in datasets:
#     file = "../final/SUPER/ALL_{}_SUPER.afreq".format(dataset)
#     print(file)
#     data[str(dataset)] = pd.read_csv(file, delimiter='\t')

for dataset in datasets:
    data[dataset] = pd.read_excel("../final/Supplementary_Data.xlsx", sheet_name=dataset)

# %%

def graph(pops, dataset, name):
    fig, ax = plt.subplots(len(pops), figsize=(25,10*len(pops)), sharex=True, sharey=True)
    plt.style.use("seaborn-darkgrid")

    for index, i in enumerate(ax):
        sns.barplot(x="Start Coordinates", y=pops[index], data=dataset, hue="ALT", dodge=False, ax=i).set(xticklabels=dataset.drop_duplicates(subset='Start Coordinates')["REF"])
        i.legend(loc='upper right')

    plt.tight_layout()
    fig.savefig('../figures/{}.png'.format(name), dpi=300)

#%%

pops = ['AFR', 'AMR', 'EAS', 'SAS', 'EUR']



for dataset in datasets:
    clinSig = pd.concat(
    [data[dataset][data[dataset][pop] >= 0.04] for pop in pops]
    , copy=False)
    display(clinSig)
    graph(pops, clinSig, "Population stratified Allele Frequency - {}".format(dataset))
# %%
