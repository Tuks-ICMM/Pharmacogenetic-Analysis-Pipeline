# %%

%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from itertools import cycle, islice
import myvariant as mv
import numpy as np
import upsetplot as pup
plt.rcParams.update({'font.size': 7})
from statannot import add_stat_annotation

# %%
# Set constants
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
pops = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
refPop = 'AFR'
compPops = ['AMR', 'EUR', 'EAS', 'SAS']
data = dict()

# %%

# Import Data
for gene in genes:
    data[gene] = pd.read_excel("../final/{}.xlsx".format(gene), sheet_name='Freq').query(" | ".join(["{} >= 0.04".format(pop) for pop in pops]))
    # data[gene] = pd.melt(data[gene], id_vars=['ID', 'POS', 'REF', 'ALT'], value_vars=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'])
    for pop in pops:
        data[gene][pop] = data[gene][pop].astype(bool)
    data[gene] = data[gene].groupby(["AFR", "AMR", "EUR", "EAS", "SAS"]).count()

# %%

for gene in genes:
    pup.plot(data[gene]['ID'], sort_by='cardinality', sort_categories_by='cardinality', show_percentages=True)
    # pup.catplot(kind='')
    plt.title("{} Intersection of alleles accros populations".format(gene), loc='left', size=10)
    # plt.savefig("../figures/{}_IntersectionPlot.jpeg".format(gene), dpi=300)

#  %%
