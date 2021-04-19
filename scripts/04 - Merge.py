# %%

# Import Dependancies
import pandas as pd
import json

# %%

# Declare Constants and Functions:
with open('../config.json') as f:
  config = json.load(f)
genes = config['locations']
clusters = config['cluster']['clusters']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
tests = ['VEP', 'Freq', 'Count', 'FishersP', 'FishersOR']

# %%

# Import data to merge:
data = dict()
for cluster in clusters:
    data[cluster] = dict()
    for gene in genes:
        data[cluster][gene] = dict()
        for test in tests:
            if cluster == 'SUB':
                if (test != 'FishersP' and test != 'FishersOR'):
                    data[cluster][gene][test] = pd.read_csv("../final/Supplementary Table/{cluster}/{gene}_{test}.csv".format(cluster=cluster, gene=gene, test=test), delimiter='\t')
                else:
                    pass
            else:
                data[cluster][gene][test] = pd.read_csv("../final/Supplementary Table/{cluster}/{gene}_{test}.csv".format(cluster=cluster, gene=gene, test=test), delimiter='\t')


# %%

# Compile a VEP + Freq sheet for analysis sake:
for cluster in clusters:
    for gene in genes:
        data[cluster][gene]['VEP_Freq'] = pd.merge(data[cluster][gene]['VEP'], data[cluster][gene]["Freq"])


# %%

# Save data to Excel
for cluster in clusters:
    for gene in genes:
        with pd.ExcelWriter("../final/{cluster}-{gene}.xlsx".format(cluster=cluster, gene=gene), engine='xlsxwriter') as writer:
            for test in data[cluster][gene]:
                df = data[cluster][gene][test]
                df.to_excel(writer, sheet_name=test, index=False, header=False)
                worksheet = writer.sheets["{}".format(test)]
                (max_row, max_col) = df.shape
                options = {
                    'columns': [{'header': column} for column in df.columns],
                    'name': test,
                    'first_column': True
                    }
                worksheet.add_table(0, 0, max_row-1, max_col-1, options)
# %%
