
#%%

import sys
import io
import json
import gzip
import requests
import time
import pandas as pd

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
params = {
    'hgvs': True,
    'CADD': True,
    'LoF': True,
    'Phenotypes': True,
    'domains': True,
    'canonical': True,
    'refseq': True
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

def chunk(seq, size):
    return [seq[pos:pos + size] for pos in range(0, len(seq), size)]

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


# %%
data = dict()
# data = read_vcf(snakemake.input["vcf"])
for location in locations:
    data[location] = read_vcf("../final/ALL_{location}.vcf.gz".format(location=location))
    # for index, row in data[location].iterrows():
    #     alleles = row['ALT'].split(",")
    #     names = row['ID'].split(";")
    #     if len(alleles) > 1 and len(names) > 1:
    #         for index, allele in enumerate(alleles):
    #             print(names)
    #             print(allele)
    #             print(index)
    #             product = row.copy()
    #             product['ALT'] = allele
    #             product['ID'] = names[index]
    #             data[location].append(product, ignore_index=True)
    #         data[location].drop(index, inplace=True)
    #     elif len(alleles) > 1 and len(names) == 1:
    #         for allele in alleles:
    #             product = row.copy()
    #             product['ALT'] = allele
    #             data[location].append(product, ignore_index=True)
    #         data[location].drop(index, inplace=True)            
    #     else:
    #         pass




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
            r = requests.post(endpoint, headers=headers, data=json.dumps(chunk), params=params)
            if r.status_code != 200:
                time.sleep(5)
            else:
                requesting = False
                decoded = r.json()
                data_received[dataset_key] = data_received[dataset_key] + decoded

# %%
# decoded = r.json()


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
frequency_data = dict()

for gene in locations:
    for pop in populations:
        supplementary[gene][pop] = 0
        frequency_data = pd.read_csv("../final/SUPER/ALL_{gene}.{pop}.afreq".format(gene=gene, pop=pop), delimiter="\t").rename(columns={'#CHROM': 'CHROM'})
        for index, row in frequency_data.iterrows():
            supplementary[gene].loc[supplementary[gene]['ID'] == row['ID'], pop] = row['ALT_FREQS']
# %%
# Rinse and repeat for clustered Frequency Analysis:
clust_freq_data = dict()

for location in locations:
    for pop in c:
        clust_freq_data[location] = dict()
        clust_freq_data[location][pop] = pd.read_csv("../final/SUPER/ALL_{location}.{pop}.afreq".format(location=location, pop=pop), delimiter="\t").rename(columns={'#CHROM': 'CHROM'})
# %%
with pd.ExcelWriter("../final/Supplementary_Data.xlsx") as writer:
    for gene in locations:
        supplementary[gene].to_excel(writer, sheet_name=gene)
    # supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)

# %%
