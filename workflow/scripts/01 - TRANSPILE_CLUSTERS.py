from os.path import join

from pandas import read_csv

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

samples = read_csv(join("input", "samples.csv"))
samples["FID"] = samples["sample_name"]
samples[["sample_name", "FID", snakemake.wildcards.cluster]].rename(
    columns={"sample_name": "IID"}
).to_csv(
    join("results", "REFERENCE", "cluster_{}.txt".format(snakemake.wildcards.cluster)),
    sep="\t",
    index=False,
)
