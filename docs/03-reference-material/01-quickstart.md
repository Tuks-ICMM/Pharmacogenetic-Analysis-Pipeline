---
title: Quickstart
layout: page
permalink: reference-material/quickstart
nav_order: 1
parent: Reference material
---

# Quickstart
{: .no_toc }

A quickreference summary of how to obtain a copy of this workflow and prep an environment for analysis.
{: .fs-6 .fw-300 }

quickstart
{: .label }

<details markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---
## Prepare working environment:
1. First you will need to download a copy of the pipeline to a location where you can configure and execute it. Navigate to our GitHub repository and retrieve the latest tag information.
2. Next, you can use `GIT` to clone a copy of the pipeline to your working environment:
    ```bash
    git clone https://github.com/Tuks-ICMM/Pharmacogenetic-Analysis-Pipeline .
    ```

    {: .normal }
    > Tags are available on our GitHub repository under the [releases](https://github.com/Tuks-ICMM/Pharmacogenetic-Analysis-Pipeline/releases) page.

## Prepare data and Metadata
1. In order to execute the _{{ site.title }}_, you will need to configure the pipeline as well as provide information about the analysis you wish to perform. This involves the following files:
    
    **Configuration files:**
      - `config/config.json` ([General configuration](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/configuration#setting-global-configuration))
    **Metadata files:**
      - `input/datasets.csv` ([Dataset declarations](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/data-requirements#datasets--dataset-files))
      - `input/samples.csv` ([Sample metadata](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/data-requirements#samples))
      - `input/locations.csv` ([Genomic location metadata](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/data-requirements#genomic-locations))
      - `input/transcripts.csv` ([Transcript selection](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/data-requirements#samples))
2. Following configuration, you will need to provide the input data files themselves.
    - `.vcf.gz` files can be compressed but must be accompanied by a tabix index file ([Discussion here](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/data-requirements#compression-and-indexing))
    - `.fasta.gz` files for reference sequences must be accompanied by a sequence dictionary file (`.dict`), a fasta index file (`.fa.gz.fai` or `fasta.gz.fai`) and a BGZIP-index (`.fa.gz.gzi`) ([Discussion here](https://tuks-icmm.github.io/Pharmacogenetic-Analysis-Pipeline/overview/configuration#reference-genomes)).

## Execute analysis
1. To execute the analysis, we need to compile our metadata and auto-generate a suitable queue-able script for the batch scheduler. To do this, you can use the `run.py` script which generates and queues a hidden generated script `.run.sh` written for your environment. For example:
    ```bash
    module load python-3.8.2
    run.py
    watch -t -d qstat
    ```

{: .normal }
> There is currently a known bug with `run.py` on some systems where the python instance does not have the requisite permissions and configuration to execute `BASH` commands. As a result, when executing the script, it will compile `.run.sh` successfully and hang without executing this script. IN such cases, the user can terminate the script (`CTRL + C`) and manually queue the generated script. In such cases, you can queue the script as follows:
>```bash
>   module load python-3.8.2
>   run.py   # CTRL + C after confirming .run.sh has been created
>   qsub .run.sh
>   watch -t -d qstat
>```


{: .normal-title }
> Support for non-PBS-Torque schedulers
>
> Currently, only support for PBS-Torque has been included in the _{{ site.title }}_. Dedicated integration of additional environments is intended for future releases. Generic snakemake profiles are currently available should you wish to manually integrate a different environments profile. Please be advised that we have made alterations to the standard Snakemake rule format to accommodate resource declarations. For more information, please see the section on the [PBS/Torque bach scheduler](/workflow-configuration/environment/pbs-torque).
>
> Documentation is available for [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and a repository of standardized profiles is available on [GitHub](https://github.com/snakemake-profiles/doc).