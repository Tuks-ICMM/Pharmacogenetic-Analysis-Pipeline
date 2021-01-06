
#%%

import sys
import io
import json
import gzip
import requests
import time
import pandas as pd

#%%

# Set constants and functions to be used:
endpoint = "https://rest.ensembl.org/vep/homo_sapiens/region/"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
params = {
    'hgvs': True,
    'CADD': True,
    'LoF': True,
    'Phenotypes': True,
    'domains': True,
    'canonical': True,
    'refseq': True,
    'transcript_id': snakemake.params.transcript_id
}

def read_vcf(path):
    """Function to import a VCF file as a Pandas dataframe. Does not include Genotype columns. CREDIT: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744

    Args:
        path (str): Path to the VCF file in question

    Returns:
        Dataframe: A Pandas Dataframe of the VCF files content.
    """
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def generate_notation(row):
    if (len(row['REF']) > len(row['ALT'])) | (len(row['REF']) == len(row['ALT'])):
        stop_coordinates = int(row['POS'])+(len(row['REF'])-1)
        return "{S}:{P}-{P2}:1/{I}".format(S=row['CHROM'], P=row["POS"], P2=stop_coordinates, I=row['ALT'])
    elif len(row['REF']) < len(row['ALT']):
        stop_coordinates = int(row['POS'])+len(row['REF'])
        return "{S}:{P}-{P2}:1/{I}".format(S=row['CHROM'], P=row["POS"], P2=stop_coordinates, I=row['ALT'])

def chunk(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def update_set(set, term, search, add_or_append=True):
    if add_or_append:
        if term in search:
            set.add(search[term])
        else:
            set.add("-")
    else:
        if term in search:
            set.append(search[term])
        else:
            set.append("-")


#%%

# data = read_vcf(snakemake.input["vcf"])
data = read_vcf(snamkemake.input['vcf'])

#%%
divided_data = chunk(data, 1000)


#%%

data_to_send = dict()
data_to_send['variants'] = list()
for chunk in divided_data:
    for index, row in chunk.iterrows():
        data_to_send['variants'].append(generate_notation(row))


#%%

requesting = True
while requesting:
    r = requests.post(endpoint, headers=headers, data=json.dumps(data_to_send), params=params)
    if r.status_code != 200:
        time.sleep(5)
    else:
        requesting=False


#%%
decoded = r.json()


#%%

i = 1

supplementary = data[['ID', 'POS', 'REF', 'ALT']]
new_columns = ["Co-Located Variant", "Transcript ID", "Transcript Strand" 'Existing Variation', 'Consequence', 'Known Alleles', "Diseases", "Biotype", "CADD_PHRED"]
for column in new_columns:
    supplementary[column] = None

for variant in decoded:
    row = supplementary.loc[supplementary['POS'] == int(variant['start'])]
    co_variants = list()
    if 'colocated_variants' in variant:
        row["Co-Located Variant"] = True
        for variant in variant['colocated_variants']:
            co_variants.append(variant['id'])
    else:
        row["Co-Located Variant"] = False
        co_variants.append("-")
    row['Existing Variation'] = ", ".join(co_variants)

    if 'transcript_consequences' in variant:
        for consequence in variant['transcript_consequences']:
            if 'canonical' in consequence:
                print(consequence['transcript_id'])
                i = i +1
                # Add fields:
                row['Consequence'] = ", ".join(consequence['consequence_terms'])
                row['Transcript ID'] = consequence['transcript_id'] if 'transcript_id' in consequence else '-'
                row["Biotype"] = consequence['biotype'] if 'biotype' in consequence else '-'
                row['CADD_PHRED'] = consequence['cadd_phred'] if 'cadd_phred' in consequence else '-'
                row['Transcript Strand'] = consequence['strand'] if 'strand' in consequence else '-'
                
                # Add phenotypes as a list:
                phenotype = set()
                if 'phenotypes' in consequence:
                    for instance in consequence['phenotypes']:
                        phenotype.add(instance['phenotype'])
            else:
                pass
# %%

supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)
