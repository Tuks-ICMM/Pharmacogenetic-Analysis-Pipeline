# %%

# Import dependancies
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import json

# %%

# Declare Constants and Functions
with open('../config.json') as f:
  config = json.load(f)
genes = config['locations']
clusters = config['cluster']['clusters']
populations = pd.read_excel("../Clusters.xlsx")
refPop = 'AFR'
compPop = ['AMR', 'EUR', 'EAS', 'SAS']


# %%
