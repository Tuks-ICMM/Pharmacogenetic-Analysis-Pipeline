# %%

%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot as pup
plt.rcParams.update({'font.size': 7})
import json

# %%

# Declare Constants and Functions
with open('../config.json') as f:
  config = json.load(f)
genes = config['locations']
clusters = config['cluster']['clusters']
populations = pd.read_excel("../Clusters.xlsx")



# %%

# Import Data:
data = dict()
for gene in genes:
    data[gene] = pd.read_excel("../final/SUPER-{}.xlsx".format(gene), sheet_name='Freq').query(" | ".join(["{} >= 0.04".format(pop) for pop in populations['SUPER'].unique()]))
    # data[gene] = pd.melt(data[gene], id_vars=['ID', 'POS', 'REF', 'ALT'], value_vars=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'])
    for pop in populations['SUPER'].unique():
        data[gene][pop] = data[gene][pop].astype(bool)
    data[gene] = data[gene].groupby(["AFR", "AMR", "EUR", "EAS", "SAS"]).count()

# %%

for gene in genes:
    pup.plot(data[gene]['ID'], sort_by='cardinality', sort_categories_by='cardinality', show_percentages=True)
    # pup.catplot(kind='')
    plt.title("{} Intersection of alleles accros populations".format(gene), loc='left', size=10)
    plt.savefig("../figures/{}_IntersectionPlot.jpeg".format(gene), dpi=300)

#  %%