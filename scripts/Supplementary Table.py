#%%

import pandas as pd

#%%

data = dict()

availData = ['CYP2A6', 'CYP2B6', 'UGT2B7']

for dataset in availData:
    data[dataset] = pd.read_csv("../data/VEP Final/{}.txt".format(dataset),  delimiter="\t")
# %%

freq = dict()

for datset in availData:
    freq[dataset] = pd.read_fwf("../data/Final/SUPER/ALL_{}_SUPER.frq.strat".format(dataset))[['SNP', 'A1', 'A2', 'CLST', 'MAF', "MAC", "NCHROBS"]]
    freq[dataset] = freq[dataset].pivot(index=["SNP", "A1", "A2"], columns="CLST", values=["MAF", "MAC", "NCHROBS"])
    #

# %%

fishers = dict()


for counter, dataset in enumerate(availData, 1):
    fishers[dataset] = pd.read_csv("../data/Fishers/Fishers_{}.csv".format(counter), delimiter="\t")
# %%
