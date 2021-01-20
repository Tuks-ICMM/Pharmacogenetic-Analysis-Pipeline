# %%

import pandas as pd

# %%

# Declare Constants and Functions:
genes = ['CYP2A6','CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
tests = ['VEP', 'Freq', 'Fishers']

# %%
# Import data to merge:
data = dict()
for gene in genes:
     data[gene] = dict()
     for test in tests:
         data[gene][test] = pd.read_csv("../final/{gene}_{test}.csv".format(gene=gene, test=test), delimiter='\t')
    

# %%
for gene in data:
    with pd.ExcelWriter("../final/{}.xlsx".format(gene)) as writer:
        for test in data[gene]:
            df = data[gene][test]
            df.to_excel(writer, sheet_name=test, index=False, header=False)

            worksheet = writer.sheets["{}".format(test)]
            (max_row, max_col) = df.shape
            column_settings = []
            for header in df.columns:
                column_settings.append({'header': header})
            worksheet.add_table(0, 0, max_row-1, max_col-1, {'columns': column_settings})
# %%
