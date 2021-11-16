import pandas as pd
from os.path import join
from peppy import Project

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Sarah Turner",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

cluster = Project(join("config", "pep.yaml")).sample_table.reset_index(drop=True)
cluster["FID"] = cluster["sample_name"]
cluster[["sample_name", "FID", snakemake.wildcards.cluster]].rename(
    columns={"sample_name": "IID"}
).to_csv(
    join("results", "REFERENCE", "cluster_{}.txt".format(snakemake.wildcards.cluster)),
    sep="\t",
    index=False,
)
