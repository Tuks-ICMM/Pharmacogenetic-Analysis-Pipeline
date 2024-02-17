from os.path import join

from pandas import read_csv

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

samples = read_csv(join("input", "samples.csv"))
samples[["sample_name", snakemake.wildcards.cluster]].rename(
    columns={"sample_name": "#IID"}
).to_csv(
    snakemake.output.metadata,
    sep="\t",
    index=False,
)
