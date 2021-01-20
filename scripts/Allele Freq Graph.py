#%% 

%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import myvariant as mv
import numpy as np
import pyupset as pup
plt.rcParams.update({'font.size': 7})

#%%
# Define constants and functions
def graph(pops, dataset, name):
    fig, ax = plt.subplots(len(pops), figsize=(25,10*len(pops)), sharex=True, sharey=True)
    plt.style.use("seaborn-darkgrid")

    for index, i in enumerate(ax):
        sns.barplot(x="POS", y=pops[index], data=dataset, hue="ALT", dodge=False, ax=i).set(xticklabels=dataset.drop_duplicates(subset='POS')["REF"])
        # sns.lineplot()
        i.legend(loc='upper right')
        i.axhline(0.04, linewidth=3, color='r', label="Clinical Significance (>4%)", ls="--")

    plt.tight_layout()
    fig.savefig('../figures/{}.png'.format(name), dpi=300)

genes = ['CYP2A6','CYP2B6', 'UGT2B7']
pops = ['AFR', 'AMR', 'EAS', 'SAS', 'EUR']
refPop = 'AFR'
compPops = ['AMR', 'EAS', 'SAS', 'EUR']
data = dict()

#%%
# Import data:
for gene in genes:
    data[gene] = pd.merge(pd.read_excel("../final/{}.xlsx".format(gene), sheet_name='Freq'), pd.read_excel("../final/{}.xlsx".format(gene), sheet_name='Fishers'), on=['ID', 'POS', 'REF', 'ALT']) 


# %%
# Compile the final, filtered dataset:
finalGraph = dict()

for gene in genes:
    toMerge = list()
    finalGraph[gene] = pd.concat(
    [data[gene][
        data[gene][pop] >= 0.04
        ] for pop in pops]
    , copy=False)
    for pop in compPops:

# %%




for gene in genes:
    clinSig = pd.concat(
    [data[gene][
        data[gene][pop] >= 0.04
        ] for pop in pops]
    , copy=False)
    graph(pops, clinSig, "Population stratified Allele Frequency - {}".format(gene))
# %%
