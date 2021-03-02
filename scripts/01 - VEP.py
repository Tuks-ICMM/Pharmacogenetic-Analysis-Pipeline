
#%%
# Import dependancies
import sys
import io
import json
import gzip
import requests
import time
import pandas as pd
from typing import Generator
import os
from typing import Tuple

# %%

# Set constants and functions to be used:
# locations = snakemake.config['locations'].keys()
with open('../config.json') as f:
  config = json.load(f)

clusters = config['cluster']['clusters']
geneSummary = dict()
for cluster in config['cluster']['clusters']:
    geneSummary[cluster] = pd.read_excel("../%s" % config["cluster"]["file"])[['ID', cluster]] 

#  %%

#  Set POST Variables and Headers:
locations=config['locations'].keys()
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
    "LoF": True,
    "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred",
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
    "LoF": True,
    "dbNSFP": "SIFT4G_score,SIFT4G_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred",
    'transcript_id': config['locations'][key]["GRCh38"]["transcript_id"]
    }
    # if key == 'CYP2A6':
    #     return params | dict(transcrcipt_id="NM_000762.6")
    # if key == 'CYP2B6':
    #     return params | dict(transcript_id="NM_000767.5")
    # if key == 'UGT2B7':
    #     return params | dict(transcript_id="NM_001074.4")

condelTests = ['SIFT', "PolyPhen"]
cutoff = {
    'SIFT': 0.15,
    'PolyPhen': 0.28,
    'Condel': 0.46
}
maximums = {
    'SIFT': 1,
    'PolyPhen': 1
}
probabilities = {
    "SIFT": pd.read_csv("../binaries/CONDEL/sift.data", delimiter='\t', names=['Score', 'Deleterious', 'Normal']),
    "PolyPhen": pd.read_csv("../binaries/CONDEL/polyphen.data", delimiter='\t', names=['Score', 'Deleterious', 'Normal'])
}


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

def generate_notation(row: pd.Series, gene: str) -> str:
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
            notation.append("{S}:{P}-{P2}:{ST}/{I}".format(S=row['CHROM'], P=row["POS"], P2=stop_coordinates, I=allele, ST=config['locations'][gene]['GRCh38']['strand']))
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

def merge(value: set) -> str:
    """Generate row value to save in df

    Args:
        value (set): The set you wish to convert

    Returns:
        str: The converted value, ready to store.
    """
    r = list(value)
    if len(r) == 1:
        return str(r[0])
    else:
        return " | ".join(value)

def CONDEL(sift: int, polyphen: int) -> Tuple[int, str]:
    """Function to calculate a weighted average score accross SIFT and PolyPhenV2 scores.

    Args:
        sift (int): A number indicating the likelyhood a variant affects protein function (Based on sequence homology and protein function)
        polyphen (int): A number indicating impact of a variant on protein structure using straightforward physical comparative methods

    Returns:
        Tuple[int, str]: A weighted score which favours variants with scores further away from each set cutoff (i.e. less ambiguity regarding the prediction) as well as a string, indicating the cutoff verdict
    """
        
    score = int()

    if sift <= cutoff['SIFT']:
        score += (1 - sift/maximums['SIFT'])*(1-probabilities['SIFT'].loc[probabilities['SIFT']['Score'] == sift, 'Normal'].values[0])
    else:
        score += (1 - sift/maximums['SIFT'])*(1-probabilities['SIFT'].loc[probabilities['SIFT']['Score'] == sift, 'Deleterious'].values[0])

    if polyphen >= cutoff['PolyPhen']:
        score += (polyphen/maximums['PolyPhen'])*(1-probabilities['PolyPhen'].loc[probabilities['PolyPhen']['Score'] == polyphen, 'Normal'].values[0])
    else:
        score += (polyphen/maximums['PolyPhen'])*(1-probabilities['PolyPhen'].loc[probabilities['PolyPhen']['Score'] == polyphen, 'Deleterious'].values[0])

    # This acocunts for number of VEP teests used. Since we are only using 2, we use 2...
    score = score/2

    print(score)

    if score >= 0.469:
        pred="deleterious"
    elif score >= 0 and score < 0.469:
        pred = "neutral"
    print (score, pred)
    
    return (score, pred)





# %%

#  Import Data:
data = dict()
for location in locations:
    data[location] = read_vcf("../final/ALL_{location}.vcf.gz".format(location=location))

 # %%

# Sub-Divide data:
data_generator = dict()
for dataset in data:
    data_generator[dataset] = chunk(data[dataset], 200)


# %%

# Compile and format request bodies:
data_to_send = dict()

# Iterate through each gene:
for dataset in data_generator:
    data_to_send[dataset] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in data_generator[dataset]:
        temp_list = list()

        # Iterate through each row in the chunk and add the HGVS notation to the list:
        for index, row in chunk.iterrows():
            temp_list.extend(generate_notation(row, dataset))
        data_to_send[dataset].append(dict(variants = temp_list))

# %%

# Perform API calls:
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

# Iterate through each response and compile its excel:
supplementary = dict()
for dataset_key, dataset in data_received.items():
    supplementary[dataset_key] = data[dataset_key][['ID', 'POS', 'REF', 'ALT']]
    new_columns = [
        'Co-Located Variant', 
        'Transcript ID', 
        'Transcript Strand', 
        'Existing Variation', 
        'Start Coordinates', 
        'Consequence', 
        'Diseases', 
        'Biotype', 
        'CADD_PHRED', 
        'LoFtool' ,
        'input',
        'SIFT4G_score',
        'SIFT4G_pred',
        'Polyphen_score',
        'Polyphen_pred',
        'CONDEL',
        'CONDEL_pred'
        ]
    for column in new_columns:
        supplementary[dataset_key][column] = "-"

    for key, chunk in data_received.items():
        for variant in chunk:
            row = supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] == int(variant['start'])]
            row['Start Coordinates'] = variant['start']
            row['Input'] = variant['input']
            
            co_variants = list()
            if 'colocated_variants' in variant:
                row["Co-Located Variant"] = True
                for colocated_variant in variant['colocated_variants']:
                    co_variants.append(colocated_variant['id'])
            else:
                row["Co-Located Variant"] = False
                co_variants.append("-")
            row['Existing Variation'] = "| ".join(co_variants)

            if 'transcript_consequences' in variant:
                Consequence = set()
                Transcript_ID = set()
                Biotype = set()
                CADD_PHRED = set()
                LoFtool = set()
                SIFT4G_score = set()
                SIFT4G_pred = set()
                Polyphen2_HVAR_score = set()
                Polyphen2_HVAR_pred = set()
                Transcript_Strand = set()
                phenotype = set()
                condel = set()
                condel_pred = set()
                for consequence in variant['transcript_consequences']:
                    # Add fields:
                    Consequence.update(consequence['consequence_terms'])
                    Transcript_ID.add(consequence['transcript_id'] if 'transcript_id' in consequence else '-')
                    Biotype.add(consequence['biotype'] if ('biotype' in consequence) else '-')
                    CADD_PHRED.add(consequence['cadd_phred'] if ('cadd_phred' in consequence) else '-')
                    # LoFtool.add(consequence['LoFtool'] if ('LoFtool' in consequence) else '-')
                    SIFT4G_score.add(consequence['sift4g_score'] if ('sift4g_score' in consequence) else '-')
                    SIFT4G_pred.add(consequence['sift_prediction'] if ('sift_prediction' in consequence) else '-')
                    Polyphen2_HVAR_score.add(str(consequence['polyphen2_hvar_score']) if ('polyphen2_hvar_score' in consequence) else '-')
                    Polyphen2_HVAR_pred.add(consequence['polyphen2_hvar_pred'] if ('polyphen2_hvar_pred' in consequence) else '-')
                    Transcript_Strand.add(str(consequence['strand'] if ('strand' in consequence) else '-'))

                    if 'sift4g_score' in consequence and 'polyphen2_hvar_score' in consequence:
                        print("yes")
                        s, p = CONDEL(consequence['sift4g_score'], consequence['polyphen2_hvar_score'])
                        condel.add(s)
                        condel_pred.add(p)
                    # Add phenotypes as a list:
                    if 'phenotypes' in consequence:
                        for instance in consequence['phenotypes']:
                            phenotype.add(instance['phenotype'])
                row['SIFT4G_score'] = merge(SIFT4G_score)
                row['SIFT4G_pred'] = merge(SIFT4G_pred)
                row['Polyphen_score'] = merge(Polyphen2_HVAR_score)
                row['Polyphen_pred'] = merge(Polyphen2_HVAR_pred)
                row['Diseases'] = merge(phenotype)
                row['Consequence'] = merge(Consequence)
                row['Transcript ID'] = merge(Transcript_ID)
                row['Biotype'] = merge(Biotype)
                row['CADD_PHRED'] = merge(CADD_PHRED)
                row['Transcript Strand'] = merge(Transcript_Strand)
                row['CONDEL'] = merge(condel)
                row['CONDEL_pred'] = merge(condel_pred)
                
            else:
                pass
            supplementary[dataset_key].loc[supplementary[dataset_key]['POS'] == int(variant['start'])] = row


# %%

# Save formatted results to CSV
for cluster in clusters:
    for gene in locations:
        supplementary[gene].to_csv("../final/Supplementary Table/{cluster}/{gene}_VEP.csv".format(cluster=cluster, gene=gene), sep='\t', index=False)
        # supplementary.to_excel(snakemake.output['excel'], sheet_name=snakemake.wildcards.location)


# %%
