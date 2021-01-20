
#%%

import sys
import io
import json
import gzip
import requests
import time
import pandas as pd
from typing import Generator
import os

# %%

# Set constants and functions to be used:
# locations = snakemake.config['locations'].keys()
locations=['CYP2A6', 'CYP2B6', 'UGT2B7']
populations = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']
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
    # 'transcript_id': snakemake.params.transcript_id
}
def generate_params(key):
    params = {
    'hgvs': True,
    'CADD': True,
    'LoF': True,
    'Phenotypes': True,
    'domains': True,
    'canonical': True,
    'refseq': True,
    # 'transcript_id': snakemake.params.transcript_id
    }
    if key == 'CYP2A6':
        return params | dict(transcrcipt_id="NM_000762.6")
    if key == 'CYP2B6':
        return params | dict(transcript_id="NM_000767.5")
    if key == 'UGT2B7':
        return params | dict(transcript_id="NM_001074.4")

#%%

# Define formulas for later use:
def read_vcf(path: str) -> pd.DataFrame:
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

def generate_notation(row: pd.Series) -> str:
    """Parses a row and returns the HGVS notation to query E! Ensembl.

    Args:
        row (Series): A row (generated from .iterrows() method) for which notation is needed.

    Returns:
        str: HGVS notation for the given variant row.
    """
    alleles = row['ALT'].split(",")
    notation = list(str())
    for allele in alleles:
        if (len(row['REF']) > len(allele)) | (len(row['REF']) == len(allele)):
            stop_coordinates = int(row['POS'])+(len(row['REF'])-1)
            notation.append("{S}:{P}-{P2}:1/{I}".format(S=row['CHROM'], P=row["POS"], P2=stop_coordinates, I=allele))
        elif len(row['REF']) < len(allele):
            stop_coordinates = int(row['POS'])+len(row['REF'])
            notation.append("{S}:{P}-{P2}:1/{I}".format(S=row['CHROM'], P=row["POS"], P2=stop_coordinates, I=allele))
        return notation

def chunk(dataset: pd.DataFrame, size: int) -> Generator:
    """Renders a generator to yield chunks of data from the original file.

    Args:
        dataset (pd.DataFrame): The dataset to divide into chunks.
        size (int): The maximum size of the resulting chunks.

    Returns:
        [type]: [description]

    Yields:
        Generator: [description]
    """
    return (dataset[pos:pos + size] for pos in range(0, len(dataset), size))



# %%
data = dict()
for location in locations:
    data[location] = read_vcf("../final/ALL_{location}.vcf.gz".format(location=location))




 # %%
data_generator = dict()
for dataset in data:
    data_generator[dataset] = chunk(data[dataset], 200)


# %%

data_to_send = dict()

# Iterate through each gene:
for dataset in data_generator:
    data_to_send[dataset] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in data_generator[dataset]:
        temp_list = list()

        # Iterate through each row in the chunk and add the HGVS notation to the list:
        for index, row in chunk.iterrows():
            temp_list.extend(generate_notation(row))
        data_to_send[dataset].append(dict(variants = temp_list))

# %%

data_received=dict()

for dataset_key, dataset in data_to_send.items():
    data_received[dataset_key] = list()
    for index, chunk in enumerate(dataset):
        requesting = True
        temp_list = list()
        while requesting:
            r = requests.post(endpoint, headers=headers, data=json.dumps(chunk), params=generate_params(dataset_key))
            if r.status_code != 200:
                time.sleep(5)
            else:
                requesting = False
                decoded = r.json()
                data_received[dataset_key] = data_received[dataset_key] + decoded

# %%
supplementary = dict()

# Iterate through each dataset and compile its excel:
for dataset_key, dataset in data_received.items():
    supplementary[dataset_key] = data[dataset_key][['ID', 'POS', 'REF', 'ALT']]
    new_columns = ['Co-Located Variant', 'Transcript ID', 'Transcript Strand', 'Existing Variation', 'Start Coordinates', 'Consequence', "Diseases", "Biotype", "CADD_PHRED"]
    for column in new_columns:
        supplementary[dataset_key][column] = "-"

    for key, chunk in data_received.items():
        for variant in chunk:
            row = supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] == int(variant['start'])]
            row['Start Coordinates'] = variant['start']
            
            co_variants = list()
            if 'colocated_variants' in variant:
                row["Co-Located Variant"] = True
                for colocated_variant in variant['colocated_variants']:
                    co_variants.append(colocated_variant['id'])
            else:
                row["Co-Located Variant"] = False
                co_variants.append("-")
            row['Existing Variation'] = ", ".join(co_variants)

            if 'transcript_consequences' in variant:
                for consequence in variant['transcript_consequences']:
                    # Add fields:
                    row['Consequence'] = ", ".join(consequence['consequence_terms'])
                    row['Transcript ID'] = consequence['transcript_id'] if 'transcript_id' in consequence else '-'
                    row["Biotype"] = consequence['biotype'] if ('biotype' in consequence) else '-'
                    row['CADD_PHRED'] = consequence['cadd_phred'] if ('cadd_phred' in consequence) else '-'
                    row['Transcript Strand'] = consequence['strand'] if ('strand' in consequence) else '-'
                    # Add phenotypes as a list:
                    phenotype = set()
                    if 'phenotypes' in consequence:
                        for instance in consequence['phenotypes']:
                            phenotype.add(instance['phenotype'])
            else:
                pass
            supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] == int(variant['start'])] = row


# %%

for gene in locations:
    supplementary[gene].to_csv("../final/{}_VEP.csv".format(gene), sep='\t', index=False)
    # supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)


# %%
