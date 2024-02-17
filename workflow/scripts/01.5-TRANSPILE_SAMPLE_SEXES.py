__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Fatima Barmania",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"
# %%

from os.path import join

from pandas import read_csv

# %%
samples = read_csv(snakemake.input.sample_annotations).rename(
    columns={"sample_name": "#IID"}
)
samples[["#IID", "sex"]].to_csv(
    snakemake.output.sample_sex_annotations,
    sep="\t",
    index=False,
)
