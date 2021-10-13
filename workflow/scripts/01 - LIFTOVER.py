#!/usr/bin/env python
"""A Python script designed to calculate frequencies, FishersExact and Variant Effect Predictions.
"""


import os
from snakemake import shell

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


# Define constant:
config = snakemake.config
wildcards = snakemake.wildcards
params = snakemake.params


# Define functions:
def directoryExists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not os.path.exists(path):
        os.makedirs(path)


print("Determining Liftover requirements now...")
if config["samples"][wildcards.sample]["refGenome"] != "GRCh38":
    shell(
        "echo 'Liftover required. All datasets have been mapped to {}'".format(
            config["samples"][wildcards.sample]["refGenome"]
        )
    ),
    shell("module load liftover"),
    if (
        config["samples"][wildcards.sample]["refGenome"] == "GRCh37"
        or config["samples"][wildcards.sample]["refGenome"] == "Hg19"
    ):
        shell("echo 'Lifting from GRCh37 to GRCh38.'"),
        directoryExists("results/LIFTOVER")
        shell(
            "module load plink-2; plink2 --vcf results/PREP/{wildcards.sample}.vcf.gz --set-all-var-ids @:#\$r-\$a --allow-extra-chr --new-id-max-allele-len 40 truncate --chr 1-22 --out results/LIFTOVER/{wildcards.sample}_PREP --export vcf-4.2 bgz --output-chr chr26"
        ),
        shell("sleep 60; tabix -p vcf results/LIFTOVER/{wildcards.sample}_PREP.vcf.gz"),
        shell(
            "module load picard-2.17.11; java {params.mem} -jar $PICARD LiftoverVcf I=results/LIFTOVER/{wildcards.sample}_PREP.vcf.gz O=results/LIFTOVER/{wildcards.sample}.vcf.gz C={params.chainFile} REJECT=results/LIFTOVER/{wildcards.sample}_REJECTED.vcf.gz R={params.ref}"
        ),
    # TODO: Add conditionals for other human reference genome builds
    else:
        print(
            "No liftover required. Dataset {} is already mapped to GRCh38.".format(
                wildcards.sample
            )
        ),
        shell("touch results/LIFTOVER/{wildcards.sample}_EXCLUDE.dat"),
        shell(
            "module load plink-1.9; plink --map input/{wildcards.sample}.map --ped input/{wildcards.sample}.ped --allow-extra-chr --chr 1-22 --recode vcf --keep-allele-order --exclude {params.exclusionList} --out results/LIFTOVER/{wildcards.sample}"
        ),
    # shell("bgzip results/LIFTOVER/{wildcards.sample}.vcf"),
    shell("sleep 1m; tabix -f -p vcf results/LIFTOVER/{wildcards.sample}.vcf.gz"),
    shell(
        "echo 'results/LIFTOVER/{wildcards.sample}.vcf.gz' >> results/LIFTOVER/merge.list"
    )
