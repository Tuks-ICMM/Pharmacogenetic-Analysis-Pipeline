# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import json
import re

# %%

# Declare Constants and Functions:
with open('../config.json') as f:
  config = json.load(f)
genes = config['locations']
clusters = config['cluster']['clusters']
populations = pd.read_excel("../Clusters.xlsx")
tests = ['VEP', 'Freq', 'Count', 'FishersP', 'FishersOR']



# %%

data = dict()

for cluster in clusters:
    data[cluster] = dict()
    for gene in genes:
        data[cluster][gene] = dict()
        for test in tests:
            if cluster == 'SUB':
                if (test != 'FishersP' and test != 'FishersOR'):
                    data[cluster][gene][test] = pd.read_excel("../final/{cluster}-{gene}.xlsx".format(cluster=cluster, gene=gene), sheet_name=test)
            else:
                data[cluster][gene][test] = pd.read_excel("../final/{cluster}-{gene}.xlsx".format(cluster=cluster, gene=gene), sheet_name=test)

# %%
for cluster in clusters:
    for gene in genes:
        for pop in populations[cluster].unique():
            data[cluster][gene]['Count']['{}_rc'.format(pop)] = data[cluster][gene]['Count']['{}_tc'.format(pop)] - data[cluster][gene]['Count']['{}_ac'.format(pop)]

# %%

longData = dict()
for cluster in clusters:
    longData[cluster] = dict()
    for gene in genes:
        dataToUse = data[cluster][gene]['Count'].rename(columns=lambda x: re.sub(r'([A-Z]{3})_([a-z]{2})', r'\2_\1', x))[[
            'ID',
            'POS',
            'REF',
            'ALT'
        ]+[
            '{type}_{pop}'.format(type=kind, pop=pop) for kind in ['tc', 'ac', 'rc'] for pop in populations[cluster].unique()
        ]]
        longData[cluster][gene] = pd.wide_to_long(dataToUse, stubnames=['ac', 'rc', 'tc'], i=['ID', 'POS', 'REF', 'ALT'], sep="_", j="Population", suffix=r"[A-Z]{3}").sort_values(by='Population')
        longData[cluster][gene] = longData[cluster][gene].rename(columns={
            "ac": "Alternate Allele Count",
            "rc": "Reference Allele Count",
            "tc": "Total Allele Observation Count"
        }).reset_index(level='Population')

# %%
# for cluster in clusters:
#     for gene in genes:
#         df = longData[cluster][gene].loc[(longData[gene]['Reference Allele Count'] != 0) & (longData[gene]['Alternate Allele Count'] != 0)]
#         plot = sns.displot(df, x='Alternate Allele Count', hue="Population", kind="kde")
#         plot.set(xscale="log", xlim=(1, None))
#         plot.fig.suptitle("{gene} Alternate Allele {type}".format(gene=gene, type="KDE Plot"))

# # %%
# for cluster in clusters:
#     for gene in genes:
#         df = longData[cluster][gene].loc[(longData[gene]['Reference Allele Count'] != 0) & (longData[gene]['Alternate Allele Count'] != 0)].reset_index()
#         plot = sns.FacetGrid(df, hue="Population")
#         plot.map(sns.barplot, "Population", "Alternate Allele Count", order=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'])
#         # plot.set(xscale="log")
#         # plot.fig.suptitle("{gene} Alternate Allele {type}".format(gene=gene, type="KDE Plot"))

# %%

subLongData = dict()

for pop in populations['SUPER'].unique():
    subLongData[pop] = dict()
    # print(pop)
    for gene in genes:
        # popSet = set(longData['SUPER'][gene][longData['SUPER'][gene]['Population'] == pop].index)
        popSet = set(populations[populations['SUPER'] == pop]['SUB'].unique())
        # print(len(popSet))
        subLongData[pop][gene] = longData['SUB'][gene][longData['SUB'][gene]['Population'].isin(popSet)]
        # display(d)
# %%
for cluster in clusters:
    if cluster == 'SUB':
        for gene in genes:
            for pop in subLongData.keys():
                df = subLongData[pop][gene].loc[(subLongData[pop][gene]['Reference Allele Count'] != 0) & (subLongData[pop][gene]['Alternate Allele Count'] != 0)].reset_index()
                df2 = pd.melt(df, id_vars=['ID', 'POS', 'REF', 'ALT', "Population"], value_vars=['Alternate Allele Count', 'Reference Allele Count'], var_name='Count type', value_name='Count')
                g = sns.FacetGrid(df2, hue="Population", col="Count type")
                g.map(sns.kdeplot, "Count")
                g.add_legend()
                g.fig.subplots_adjust(top=0.8)
                g.fig.suptitle("{gene} Alternate Allele {type} | {pop} {clst}".format(pop=pop, gene=gene, type="KDE Plot", clst=cluster))
                # plot = sns.displot(df, x='Alternate Allele Count', y="Reference Allele Count", hue='Population')
                # plot.fig.suptitle("{} Alternate Allele Kernel Density".format(gene))
                g.savefig("../figures/KDE/{cluster}/{pop} - {gene} Alternate Allele {type}.jpeg".format(pop=pop, cluster=cluster, gene=gene, type="KDE Plot"), dpi=300)

    else:
        for gene in genes:
            df = longData[cluster][gene].loc[(longData[cluster][gene]['Reference Allele Count'] != 0) & (longData[cluster][gene]['Alternate Allele Count'] != 0)].reset_index()
            df2 = pd.melt(df, id_vars=['ID', 'REF', 'ALT', "Population"], value_vars=['Alternate Allele Count', 'Reference Allele Count'], var_name='Count type', value_name='Count')
            g = sns.FacetGrid(df2, hue="Population", col="Count type")
            g.map(sns.kdeplot, "Count")
            g.add_legend()
            g.fig.subplots_adjust(top=0.8)
            g.fig.suptitle("{gene} Alternate Allele {type} | {clst}".format(gene=gene, type="KDE Plot", clst=cluster))
            # plot = sns.displot(df, x='Alternate Allele Count', y="Reference Allele Count", hue='Population')
            # plot.fig.suptitle("{} Alternate Allele Kernel Density".format(gene))
            g.savefig("../figures/KDE/{cluster}/{gene} Alternate Allele {type}.jpeg".format(cluster=cluster, gene=gene, type="KDE Plot"), dpi=300)


# %%