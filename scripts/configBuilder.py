# %%
import pandas as pd
import json
# %%

genes = pd.read_excel("./gene_summary.xlsx", sheet_name='Simplified Masters')


# %%

j = dict()
for index, item in genes.iterrows():
    j[item['Gene']] = dict()
    j[item['Gene']]['GRCh38'] = dict()
    j[item['Gene']]['GRCh38']['transcript_id'] = item['Ref Transcript']
    j[item['Gene']]['GRCh38']['chromosome'] = item['Chr']
    j[item['Gene']]['GRCh38']['strand'] = item['Strand']
    j[item['Gene']]['GRCh38']['to'] = item['Start']
    j[item['Gene']]['GRCh38']['from'] = item['Stop']
# %%

with open("config.location.json", "w") as write_file:
    json.dump(j, write_file)

# %%
