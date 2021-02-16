# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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

for gene in genes:
    for pop in populations:
        data[gene]['Fishers']['{}_rc'.format(pop)] = data[gene]['Fishers']['{}_tc'.format(pop)] - data[gene]['Fishers']['{}_ac'.format(pop)]


# %%
for gene in genes:
    sns.displot(data=data[gene]['Fishers'].loc[data[gene]['Fishers']['AFR_ac'] != 0], x='AFR_ac')

# %%

longData = dict()

for gene in genes:
    dataToUse = data[gene]['Fishers'].rename(columns={
        'AFR_tc': 'tc_AFR',
        'AFR_ac': 'ac_AFR',
        'AFR_rc': 'rc_AFR',
        'AMR_tc': 'tc_AMR',
        'AMR_ac': 'ac_AMR',
        'AMR_rc': 'rc_AMR',
        'EUR_tc': 'tc_EUR',
        'EUR_ac': 'ac_EUR',
        'EUR_rc': 'rc_EUR',
        'EAS_tc': 'tc_EAS',
        'EAS_ac': 'ac_EAS',
        'EAS_rc': 'rc_EAS',
        'SAS_tc': 'tc_SAS',
        'SAS_ac': 'ac_SAS',
        'SAS_rc': 'rc_SAS',
    })[[
        'ID', 
        'REF', 
        'ALT',
        'tc_AFR',
        'ac_AFR',
        'rc_AFR',
        'tc_AMR',
        'ac_AMR',
        'rc_AMR',
        'tc_EUR',
        'ac_EUR',
        'rc_EUR',
        'tc_EAS',
        'ac_EAS',
        'rc_EAS',
        'tc_SAS',
        'ac_SAS',
        'rc_SAS'
        ]]
    longData[gene] = pd.wide_to_long(dataToUse, stubnames=['ac', 'rc', 'tc'], i=['ID', 'REF', 'ALT'], sep="_", j="Population", suffix=r"[A-Z]{3}").sort_values(by='Population')
    longData[gene] = longData[gene].rename(columns={
        "ac": "Alternate Allele Count",
        "rc": "Reference Allele Count",
        "tc": "Total Allele Observation Count"
    }).reset_index(level='Population')

# %%

for gene in genes:
    df = longData[gene].loc[(longData[gene]['Reference Allele Count'] != 0) & (longData[gene]['Alternate Allele Count'] != 0)]
    plot = sns.displot(df, x='Alternate Allele Count', hue="Population", kind="kde")
    plot.set(xscale="log", xlim=(1, None))
    plot.fig.suptitle("{gene} Alternate Allele {type}".format(gene=gene, type="KDE Plot"))

# %%

for gene in genes:
    df = longData[gene].loc[(longData[gene]['Reference Allele Count'] != 0) & (longData[gene]['Alternate Allele Count'] != 0)].reset_index()
    plot = sns.FacetGrid(df, hue="Population")
    plot.map(sns.barplot, "Population", "Alternate Allele Count", order=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'])
    # plot.set(xscale="log")
    # plot.fig.suptitle("{gene} Alternate Allele {type}".format(gene=gene, type="KDE Plot"))

# %%
for gene in genes:
    df = longData[gene].loc[(longData[gene]['Reference Allele Count'] != 0) & (longData[gene]['Alternate Allele Count'] != 0)].reset_index()
    df2 = pd.melt(df, id_vars=['ID', 'REF', 'ALT', "Population"], value_vars=['Alternate Allele Count', 'Reference Allele Count'], var_name='Count type', value_name='Count')
    g = sns.FacetGrid(df2, hue="Population", col="Count type")
    g.map(sns.kdeplot, "Count")
    g.add_legend()
    # plot = sns.displot(df, x='Alternate Allele Count', y="Reference Allele Count", hue='Population')
    # plot.fig.suptitle("{} Alternate Allele Kernel Density".format(gene))

# %%
