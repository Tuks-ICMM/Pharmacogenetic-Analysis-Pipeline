# %%
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# %%


data = sns.load_dataset('titanic')
data = data[data['fare'] >= 60].reset_index()
data['index'] = data['index'].astype(str)
# %%

fig, ax = plt.subplots(len(data["class"].unique()), figsize=(25,10*len(data["class"].unique())), sharex=True, sharey=True)
plt.style.use("ggplot")

for index, i in enumerate(ax):
    sns.barplot(x='index', y='fare', hue='class', dodge=False, data=data, ax=i)
    sns.scatterplot(x='index', y='fare', data=data.query("age >= 25"), zorder=10, color='red', markers=["*"], ax=i)
    
# %%
