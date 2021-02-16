# %%

# Import Dependancies
import pandas as pd

# %%

# Declare Constants and Functions:
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
tests = ['VEP', 'Freq', 'FishersP', 'FishersOR']

# %%

# Import data to merge:
data = dict()
for gene in genes:
     data[gene] = dict()
     for test in tests:
         data[gene][test] = pd.read_csv("../final/Supplementary Table/{gene}_{test}.csv".format(gene=gene, test=test), delimiter='\t')
    

# %%

# Save data to Excel
for gene in genes:
    with pd.ExcelWriter("../final/{}.xlsx".format(gene)) as writer:
        for test in data[gene]:
            df = data[gene][test]
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
