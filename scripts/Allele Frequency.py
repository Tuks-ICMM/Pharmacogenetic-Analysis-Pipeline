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


for dataset in datasets:
    file = "../final/SUPER/ALL_{}_SUPER.afreq".format(dataset)
    print(file)
    data[str(dataset)] = pd.read_csv(file, delimiter='\t')
# %%

def graph(pops, dataset, name):
    fig, ax = plt.subplots(len(pops), figsize=(25,10*len(pops)), sharex=True, sharey=True)
    plt.style.use("seaborn-darkgrid")

    for index, i in enumerate(ax):
        sns.barplot(x="start_coord", y=pops[index], data=dataset, hue="A2", dodge=False, ax=i).set(xticklabels=data["A1"])

    plt.tight_layout()
    plt.savefig('../Figures/{}.png'.format(name), dpi=300)

#%%

pops = ['AFR', 'AMR', 'EAS', 'SAS', 'EUR']

for dataset in datasets:
    graph(pops, data, "Population stratified Allele Frequency | {}".format(dataset))
# %%
