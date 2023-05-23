---
title: Data Requirements
permalink: overview/data-requirements
layout: page
nav_order: 1
parent: Workflow
---

# Data Requirements
{: .no_toc}

What files do I need to provide as inputs for this analysis?
{: .fs-6 .fw-300 }

Datasets
{: .label }

Sample Metadata
{: .label }

Genomic Location Metadata
{: .label }

Transcript Selection
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

This page lists the informational requirements needed to execute the _{{ site.title }}_. Below is a simplified overview diagram indicating the full list of required information. For more information, you may consult the relevant section below which contains more advanced explanations, discussions and technical documentation for each requirement and any formatting involved.

## Overview

```mermaid
---
title: Data Requirements
---
flowchart TB
  subgraph Standard ["Standard Resources"]
            genomeFasta[/"Reference Genome GRCh38 (FASTA)"/]
        end
        subgraph projectSpecific ["Project specific data"]
            %% Use LR to invert axis set by parent to effectively force relative "TB"
            direction LR
            subgraph data ["Variant input data"]
                datasetFiles[/"Datasets (VCF)"/]
            end
            subgraph metadata ["Analysis metadata"]
                %% Use LR to invert axis set by parent to effectively force relative "TB"
                direction LR

                datasetMeta[/"Datasets metadata (CSV)"/]
                locationMeta[/"Genomic location metadata (CSV)"/]
                sampleMeta[/"Sample metadata (CSV)"/]
                transcriptMeta[/"Transcript metadata (CSV)"/]
            end
        end
```

## Reference Genome

Due to the nature of reference sequences, they are not included as pipeline inputs and are, as a result, not handled in this section. For additional information, please consult the Reference Genome section on the pipeline configuration page.

## Datasets & Dataset Files

Please provide all input datasets in the form of _variant-call-format_ or `.vcf` files. The latest version of the VCF specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.3.pdf).

### Compression and Indexing

Due to the nature of bioinformatics and genomics, datasets are often quite large in uncompressed form. Users are welcome to compress their data files for additional performance and administrative ease-of-use.

If you wish to compress your VCF files, please provide the following files as input:

- [x] BGZIP-compressed VCF file (`.vcf.gz` or `vcf.bgz`)
- [x] Tabix Index (`.vcf.gz.tbi` or `.vcf.bgz.tbi`)

{: .normal }

> This pipeline is designed to accept `.vcf.gz` files produced by **Block Compression (BGZIP)**. This is a non-standard type of compression which is not typically the default on Windows or MacOS. It is used to compress `.vcf` files in a series of blocks or chunks and can be done using many popular bioinformatics tools.
>
> Normally, block-compression alone would only make your data file smaller. To facilitate more efficient usage of computational resources, you can also create a **Tabix Index**. This is an accompanying index file to BGZIP-compressed `.vcf.gz` files which contains an index indicating the bounds of each compression block relative to the genomic coordinates and variant IDs in the dataset, making targeted decompression much more efficient.
>
> Both block-compression and tabix indexing are provided as part of [SamTools](http://www.htslib.org/doc/bgzip.html).

## Metadata Declarations

To run the _{{ site.title }}_, you will need to provide some additional contextual information. All metadata is provided in the form of appropriately named ` .csv` files located in the `input` directory.

{: .normal-title }

> Case sensitivity
>
> The following metadata declaration files use _**case-sensitive column names**_.

---
### Datasets

The `datasets.csv` file allows you to declare datasets and provide the necessary dataset-level information for use in this pipeline.

#### Data requirements

<dl class="def-wide">
  <dt>dataset_name <code>&lt;str&gt;</code></dt>
  <dd>The name of the dataset. This value will be used as a universal accessor for that dataset and any information relating to it. This means that any output files will use this value to determine things like filenames, etc. It is also used to connect other metadata to this dataset computationally, E.g. sample-level information.
  
  <br><strong><i>E.g. <code>1000G</code></i></strong></dd>
  
  <dt>reference_genome <code>&lt;str&gt;</code></dt>
  <dd>An <code>enum</code> indicating which reference genome version this dataset has been called on.
  
  <br><strong><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></strong></dd>
  
  <dt>file <code>&lt;file_path&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis.
  
  <br><strong><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></strong></dd>
</dl>

#### `datasets.csv` data example

| **dataset_name** | **reference_genome** | **file**                                                    |
| :--------------- | :------------------- | :---------------------------------------------------------- |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |

---
### Samples

The `samples.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

#### Data requirements

<dl class="def-wide">
  <dt>sample_name <code>&lt;str&gt;</code></dt>
  <dd>The ID of the sample. this should correspond to the sample ID's provided in the provided <code>.vcf</code> file. 
  
  <br><strong><i>E.g. <code>HG002</code></i></strong></dd>
  
  <dt>dataset <code>&lt;enum [dataset_name]&gt;</code></dt>
  <dd>The name of the dataset this sample belongs to. This value should correspond to the provided dataset ID listed in <code>datasets.csv</code> 
  
  <br><strong><i>E.g. <code>1000g</code></i></strong></dd>
  
  <dt><code>* &lt;str&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis.
  
  <br><strong><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></strong></dd>
</dl>

#### `samples.csv` data example

| **sample_name** | **dataset** | **SUPER** | **SUB** |
| :-------------- | :---------- | :-------- | :------ |
| HG002           | HG002       | `EUR`     | `GBR`   |
| HG002           | HG003       | `AFR`     | `GWD`   |
| HG002           | HG004       | `SAS`     | `GIH`   |

---
### Genomic Locations

The `locations.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

#### Data requirements

<dl class="def-wide">
  <dt>location_name <code>&lt;str&gt;</code></dt>
  <dd>The ID of a gene or, if not a studied gene region, a unique identifier to reference this genomic coordinate window.
  
  <br><strong><i>E.g. <code>CYP2A6</code></i></strong></dd>
  
  <dt>chromosome <code>&lt;enum &lt;int [0-24]&gt; &gt;</code></dt>
  <dd>The chromosome number on which the above genomic region can be found.
  
  <br><strong><i>E.g. <code>19</code></i></strong></dd>

  <dt>start <code>&lt;int&gt;</code></dt>
  <dd>The start coordinates for the genomic window.
  
  <br><strong><i>E.g. <code>40842850</code></i></strong></dd>
  
  <dt>stop <code>&lt;int&gt;</code></dt>
  <dd>The stop coordinates for the genomic window.
  
  <br><strong><i>E.g. <code>1000g</code></i></strong></dd>
  
  <dt>strand <code>&lt;enum [-1,1]&gt;</code></dt>
  <dd>The strand on which the genomic region can be found, where <code>1</code> denotes the forward strand and <code>-1</code> denotes the reverse strand.
  
  <br><strong><i>E.g. <code>-1</code></i></strong></dd>
</dl>

#### `locations.csv` data example

| **location_name** | **chromosome** | **start** | **stop**  | **strand** |
| :---------------- | :------------- | :-------- | :-------- | :--------- |
| CYP2A6            | 19             | 40842850  | 40851138  | -1         |
| CYP2B6            | 19             | 40988570  | 41021110  | 1          |
| UGT2B7            | 4              | 69045214  | 69112987  | 1          |

---
### Transcripts

The `transcripts.csv` file allows you to declare which transcripts you would like to use when performing variant-effect-prediction.

During the execution of the _{{ site.title }}_, variant-effect-prediction (VEP) is performed using a publicly accessible VEP query API by E! Ensembl. Currently, the API returns multiple VEP predictions based on any transcripts that are present at a given genomic location. Users are able to provide a <code>transcripts.csv</code> input file to declare a list of transcripts per genomic-region they would like to consider for this analysis. 

{: .normal-title }
> Transcript IDs
>
>Please use transcripts listed on the [E! Ensembl Database](https://www.ensembl.org/)

{: .normal-title }
> Multiple Transcripts
>
> If more than one transcript is provided for a given genomic region, we will attempt to match the transcripts available in the order that is provided from top to bottom. The first successful VEP transcript match between the users selection and that provided by E! Ensembl will be selected, and if no transcripts provided are available, the first available transcript result will be selected.

#### Data requirements

<dl class="def-wide">
  <dt>gene_name <code>&lt;enum [str]&gt;</code></dt>
  <dd>The name of the gene a transcript describes. This key should match the gene or region name provided in the <code>locations.csv</code> file. 
  
  <br><strong><i>E.g. <code>HG002</code></i></strong></dd>
  
  <dt>transcript_id <code>&lt;str&gt;</code></dt>
  <dd>The name of the transcript in question. This value will be used to query the E! Ensembl database when performing variant-effect-prediction. 
  
  <br><strong><i>E.g. <code>NM_000762.6</code></i></strong></dd>
</dl>

#### `transcripts.csv` data example

| **gene_name** | **transcript_id**   |
| :------------ | :------------------ |
| CYP2A6        | NM_000762.6         |
| CYP2A6        | ENST00000600495.1   |
| CYP2A6        | ENST00000596719.5   |
| CYP2A6        | ENST00000599960.1   |
| CYP2B6        | NM_000767.5         |
| CYP2B6        | ENST00000593831.1   |
| CYP2B6        | ENST00000598834.2   |
| CYP2B6        | ENST00000597612.1   |
| CYP2B6        | ENST00000594187.1   |
| UGT2B7        | NM_001074.4         |
| UGT2B7        | ENST00000508661.5   |
| UGT2B7        | ENST00000622664.1   |
| UGT2B7        | ENST00000502942.5   |
| UGT2B7        | ENST00000509763.1   |
