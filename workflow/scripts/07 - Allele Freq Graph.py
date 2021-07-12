#%% 

%matplotlib inline

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import json
# plt.rcParams.update({'font.size': 25})
sns.set(font_scale=3)

#%%

# Declare Constants and Functions:
with open('../config.json') as f:
  config = json.load(f)
genes = config['locations']
pops = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
refPop = 'AFR'
compPops = ['AMR', 'EUR', 'EAS', 'SAS']

# %%

# Define constants and functions
def graph(pops, dataset, markers, name):
    fig, ax = plt.subplots(len(pops), figsize=(25,10*len(pops)), sharex=True, sharey=True)
    plt.style.use("ggplot")
    # data = pd.DataFrame(dataset['ALT'].unique(), columns=['ALT'])
    # data['color'] = pd.Series(plt.get_cmap('viridis', len(data['ALT'])).colors)
    # display(data)
    # fig.suptitle(name)

    for index, i in enumerate(ax):

        g = sns.barplot(x='ID', y=pops[index], data=dataset, ax=i, hue='ALT', dodge=False)
        # for indx, row in markers.iterrows():
            # print(row['POS'])
            #g.text(row['POS'], row[pops[index]], "*", color='red')
        # sns.scatterplot(x='POS', y=pops[index], data=markers, zorder=10, color='red', markers=["*"], ax=i)
        # for indx, item in enumerate(dataset['ALT'].unique()):
        #     data = dataset[dataset['ALT'] == item]
        # i.bar(data['POS'], data[pops[index]], align='center', color=(lambda alt:  data.loc[data['ALT'] == alt, 'color']), label=item)#.set(xticklabels=dataset.drop_duplicates(subset='POS')["REF"])
        # i.plot(markers['POS'], markers[pops[index]], marker="*", color='r', linestyle='None', markersize=20.0)
        # i.set_yscale("log")
        i.tick_params(axis='x', which='major', labelrotation = 90)
        i.tick_params(axis='y', which='major', labelrotation = 45, fontsize=15)
        # i.set_xticklabels(dataset['ID'], rotation = 45, ha='right')
        i.legend(loc='upper right', fontsize=25, title_fontsize=30, framealpha=0.5)
        i.axhline(0.04, linewidth=5, color='r', label="Clinical Significance (>=4%)", ls="--")
        i.axhline(0.01, linewidth=2, color='gray', label="Allele Cuttoff (>=1%)", ls="--")

    plt.tight_layout()
    # fig.savefig('../figures/{}.png'.format(name), dpi=300)

#%%
# Import data:
data = dict()

for gene in genes:
    data[gene] = pd.merge(pd.read_excel("../final/SUPER-{}.xlsx".format(gene), sheet_name='Freq'), pd.read_excel("../final/SUPER-{}.xlsx".format(gene), sheet_name='FishersP'), on=['ID', 'POS', 'REF', 'ALT'])
    data[gene]['POS'] = data[gene]['POS'].astype(str)

# %%

# Filter data to Fishers Significance:
clinSig = dict()
clinSigP = dict()

for gene in genes:
    clinSig[gene] = data[gene].query(" | ".join(["{} >= 0.04".format(pop) for pop in compPops])).reset_index().drop('index', axis='columns')
    
    # pd.concat(
    # [data[gene][
    #     data[gene][pop] >= 0.04
    #     ] for pop in pops]
    # , copy=False, join='inner')
    # for pop in compPops:
    
    clinSigP[gene] = clinSig[gene].query(" | ".join(["AFR_T_{pop} >= 0.05".format(pop=pop) for pop in compPops]))
    clinSigP[gene] = clinSig[gene].query(" | ".join(["AFR_G_{pop} >= 0.05".format(pop=pop) for pop in compPops]))
    clinSigP[gene] = clinSig[gene].query(" | ".join(["AFR_L_{pop} >= 0.05".format(pop=pop) for pop in compPops]))
    # graph(pops, clinSig[gene].set_index("ID"), clinSigP[gene].set_index("ID"), "Population stratified Allele Frequency - {}".format(gene))
  # %%

  # Graph each gene:
for gene in genes:
    graph(pops, clinSig[gene].sort_values(by='POS'), clinSigP[gene].sort_values(by='POS'), "Population stratified Allele Frequency - {}".format(gene))
  # %%
