# SnakeMake Pharmacogenetics pipeline

A pipeline to perform frequency analysis on a dataset of African populations. The dataset comprises of individual genotype data obtained from 1000 Genomes Project (1000g), African Genome Variation Project (AGVP), South African Human Genome Project (SAHGP). This is the initial steps towards fully automating the pipeline to accomodate any datasets.

## Requirements:
This has been made using a python-based DSL called [Snakemake](https://snakemake.readthedocs.io/en/stable/)and coded to run on a PBS/Torque environment using the `qsub` command (Hard-coded into each process). As such, it needs to be run on a server with the appropriate binaries and batch scheduling software. Below is a list of software used by this pipeline:
- PBS/Torque batch scheduler
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [PLINK-1.9](https://www.cog-genomics.org/plink2)
- [VCF-Tools](https://vcftools.github.io/index.html)
- [liftOverPlink](https://github.com/sritchie73/liftOverPlink)(Binaries contained within this repo. _**Update at own risk!**_)
  - [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)(Required dependancy for liftOverPlink)
- [e! Ensembl VEP API](https://www.ensembl.org/info/docs/tools/vep/index.html)

## Instructions:
1. use the `cd` command to navigate to the root repository directory containing the `Snakefile`. (This file contains all of the instructions for each step of the analysis)
2. To start the pipeline and produce the default list of files, simply call `snakemake` on the command line.
3. To generate a runtime report, detailing figures produced and performance-related numbers, use the `snakemake --report` command. The HTML file produced is completely self-contained and can be shared as needed. You can view it using any web browser such as firefox or Google Chrome, etc.

## ToDo:
> Ordered by priority:
- [ ] Segment the processes within iterative loops to handle files passed iteratively and dynamically no matter how many are passed through numbers-wise.
- [X] ~~Wire in command-line argument support to make pipeline sharable.~~ _(Made redundant by Snakemake implementation)_
- [ ] Integrate variant-effect-prediction tools in-house as final pipeline portion ([SIFT](https://sift.bii.a-star.edu.sg/), [PolyPhen](http://genetics.bwh.harvard.edu/pph2/), [CONDEL] (https://omictools.com/condel-tool), etc)
- [ ] **(OPTIONAL)** If time permits, try and integrate graph generation into pipeline using Python and [matplotlib](https://matplotlib.org/).