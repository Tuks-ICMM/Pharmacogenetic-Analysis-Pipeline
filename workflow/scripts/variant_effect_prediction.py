#!/usr/bin/env python
"""
A Python script designed to run Variant Effect Prediction calculations
and API calls as well as save the results to file.

[NEEDS] WORKFLOW_RUNTIME

"""
# %%
############ IMPORT DEPENDANCIES ############
#############################################

from json import dumps
from os.path import join
from time import sleep

import pandas as pd
from common.common import (
    chunk,
    directory_exists,
    generate_notation,
    generate_params,
    merge,
    read_vcf,
    save_or_append_to_excel,
)
from common.condel_score import condel_weighted_score
from requests import post

# from workflow.scripts.entities.VariantConsequence import VariantConsequenceResult

# [SET] Pandas Chained-Assignment off to prevent writing changes to transient DF copies:
pd.set_option("chained_assignment", None)


__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

# %%
############ IMPORT DATA ############
#####################################

# [IMPORT] sample-level data:
SAMPLE_DATA = pd.read_csv(join("..", "..", "input", "samples.csv"))

# [IMPORT] gene-level data:
GENE_DATA = pd.read_csv(join("..", "..", "input", "locations.csv"))

# [IMPORT] transcript selections:
TRANSCRIPT_DATA = pd.read_csv(join("..", "..", "input", "transcripts.csv"))


# %%
############ SET CONSTANT VARIABLES ############
################################################

# [SET] cluster names:
CLUSTERS = set(SAMPLE_DATA.keys())
CLUSTERS.remove("sample_name")
CLUSTERS.remove("dataset")

# [SET] gene names:
GENES = GENE_DATA["location_name"].unique().tolist()

# [SET] population metadata:
REFERENCE_POPULATION = "AFR"
COMPARISON_POPULATIONS = ["AMR", "EUR", "EAS", "SAS"]
ALL_POPULATIONS = ["AFR", "AMR", "EUR", "EAS", "SAS"]

# [SET] variant effect prediction metadata
ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region/"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# [SET] API payload categories for results:
VARIANT_TYPE_CATEGORIES = [
    "colocated_variants",
    "transcript_consequences",
    "regulatory_feature_consequences",
]
# %%
############ READ IN VCF FILE AS BASE ############
##################################################

RAW_DATA_TEMPLATE = dict()
for gene in GENES:
    RAW_DATA_TEMPLATE[gene] = read_vcf(
        join("..", "..", "results_OLD_BACKUP", "FINAL", "super-population", f"ALL_{gene}.vcf.gz")
    )[["CHROM", "POS", "ID", "REF", "ALT"]]
    RAW_DATA_TEMPLATE[gene]["QUERY"] = None
    # Lets create a multi-index to access specific variant records, more easily:
    RAW_DATA_TEMPLATE[gene].set_index(["CHROM", "POS", "REF", "ALT"], inplace=True)

# %%
############ SUB_DIVIDE DATA ############
#########################################
DATA_IN_CHUNKS = dict()
for gene, data in RAW_DATA_TEMPLATE.items():
    chunk_var = chunk(data, 200)
    DATA_IN_CHUNKS[gene] = chunk_var


# %%k
############ COMPILE E! ENSEMBL REQUEST PAYLOADS ############
#############################################################
# Compile and format request bodies:
API_QUERY_STATEMENTS = dict()
# Iterate through each gene:
for gene in GENES:
    # Nest, lets create a placeholder list that we can append
    # the query strings to in preparation to make the API calls
    API_QUERY_STATEMENTS[gene] = list()

    # Iterate through each n-sized chunk generated:
    for chunk in DATA_IN_CHUNKS[gene]:
        temp_list = list()

        # Iterate through each row in the chunk and add the HGVS notation to the list:
        for index, row in chunk.iterrows():
            # Here we extract the multiindex key for this spesific variant record.
            # This is done because we are going to build a second indexable column. This
            # column will be to track the E! Ensembl query string and will also serve as
            # the new index when annotating teh DataFrame with the API payloads.
            CHROM, POS, REF, ALT = row.name

            # Here, we generate the query string for the E! Ensembl API:
            ensembl_query_string: str = generate_notation(
                row,
                GENE_DATA.loc[GENE_DATA["location_name"] == gene, "strand"].item(),
            )

            # Add the generated notation to the growing temp list of variants in this chunk of data
            RAW_DATA_TEMPLATE[gene].loc[
                (CHROM, POS, REF, ALT), "QUERY"
            ] = ensembl_query_string
            temp_list.append(ensembl_query_string)

        # Add the chunks worth of formatted variant query strings into the payload object
        API_QUERY_STATEMENTS[gene].append(dict(variants=temp_list))

# %%
########### PERFORM E! ENSEMBL API CALLS ############
#####################################################
API_RESPONSE = dict()
for gene in API_QUERY_STATEMENTS.keys():
    API_RESPONSE[gene] = list()
    for index, chunk in enumerate(API_QUERY_STATEMENTS[gene]):
        REQUESTING = True
        while REQUESTING:
            r = post(
                ENDPOINT,
                headers=HEADERS,
                data=dumps(chunk),
                params=generate_params(gene, True),
            )
            if not r.ok:
                print(str(r.reason))
                print(r.json())
                sleep(2)
            else:
                decoded = r.json()
                API_RESPONSE[gene].extend(decoded)
                REQUESTING = False
# %%

import pickle

with open("ENSEMBL_API_PAYLOAD_WITH_PICK_FLAGS.pkl", "wb") as file:
    pickle.dump(API_RESPONSE, file)

# %%

############ OPEN 'CANNED' E! ENSEMBL API CALLS ############
############################################################
# import pickle

with open("ENSEMBL_API_PAYLOAD.pkl", "rb") as file:
    API_RESPONSE = pickle.load(file)

# %%
# TODO: Improve this piece of code
# DATA = dict()
# for gene, genomic_location_responses in API_RESPONSE.items():
#     DATA[gene] = RAW_DATA_TEMPLATE[gene].set_index("QUERY")
#     # [["CHROM", "ID", "POS", "REF", "ALT"]]
#     # DATA[gene].set_index(["CHROM", "POS", "REF", "ALT"])
#     new_columns = []
#     for column in new_columns:
#         DATA[gene][column] = None

#     for response in genomic_location_responses:
#         resp: VariantConsequenceResult = VariantConsequenceResult.from_dict(response)
    #     consequence = resp.transcript_consequences
    #     DATA[gene].loc[
    #         (resp.input), "colocated variant start"
    #     ] = resp.generateDataFrameRow()


# %%
############ COMPILE E! ENSEMBL API RESPONSES ############
##########################################################
# Abstract response parsing to functions
##########################################################
DATA = dict()
for gene in GENES:
    DATA[gene] = dict()
    for category in VARIANT_TYPE_CATEGORIES:
        # Filter records per category so that we can group all related variant types in groups:
        filtered_subset = list(
            filter(lambda record: category in record, API_RESPONSE[gene])
        )
        if category == "colocated_variants":
            for record in filtered_subset:
                for colocated_variant in record["colocated_variants"]:
                    if "frequencies" in colocated_variant:
                        colocated_variant["frequencies_gnomad"] = colocated_variant[
                            "frequencies"
                        ][record["allele_string"].split("/")[1]]
                        del colocated_variant["frequencies"]

        if filtered_subset:
            normalized_dataframe = pd.json_normalize(
                filtered_subset,
                category,
                record_prefix=f"{category}.",
                meta=[
                    "assembly_name",
                    "seq_region_name",
                    "most_severe_consequence",
                    "id",
                    "start",
                    "strand",
                    "variant_class",
                    "input",
                    "end",
                ],
            )
            DATA[gene][category] = pd.merge(
                RAW_DATA_TEMPLATE[gene].reset_index(),
                normalized_dataframe,
                how="left",
                left_on="QUERY",
                right_on="input",
            )

            # DATA[gene][category] = merged_dataframe.reindex(
            #     [
            #         "CHROM",
            #         "POS",
            #         "REF",
            #         "ALT",
            #         "assembly_name",
            #         "seq_region_name",
            #         "most_severe_consequence",
            #         "input",
            #         "id",
            #         "start",
            #         "end",
            #         "strand",
            #         "variant_class",
            #     ], axis='columns'
            # )

# %%


# %%
############ SAVE TO EXCEL ############
#######################################

# [SAVE] VEP analysis to excel
for cluster in CLUSTERS:
    directory_exists(join("..", "..", "results", "FINAL"))
    for gene in GENES:
        for category in DATA[gene].keys():
            save_or_append_to_excel(
                DATA[gene][category], cluster, gene, category, "replace"
            )
# %%
