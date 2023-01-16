---
title: Data Requirements
permalink: overview/data-requirements
nav_order: 3
layout: page
---

# Data Requirements

This page lists the informational requirements needed to execute the _{{ site.title }}_.

## Reference Genome

## Datasets

### Data files

## Metadata Declarations

To run the _{{ site.title }}_, you will need to provide some additional contextual information. All metadata is provided in the form of appropriately named ` .csv`` files located in the  `input` directory.

{: .note-title }

> Case sensitivity
>
> The following metadata declaration files use _**case-sensitive column names**_.

### Datasets

The `datasets.csv` file allows you to declare datasets and provide the necessary dataset-level information for use in this pipeline.

#### Data requirements

<dl>
  <dt>dataset_name <code>&lt;str&gt;</code></dt>
  <dd>The name of the dataset. This value will be used as a universal accessor for that dataset and any information relating to it. This means that any output files will use this value to determine things like filenames, etc. It is also used to connect other metadata to this dataset computationally, E.g. sample-level information.
  
  E.g. <code>1000G</code></dd>
  
  <dt>reference_genome <code>&lt;str&gt;</code></dt>
  <dd>An <code>enum</code> indicating which reference genome version this dataset has been called on.
  
  E.g. <code>GRCh37</code> or <code>GRCh38</code></dd>
  
  <dt>file <code>&lt;file_path&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis.
  
  E.g. <code>GRCh37</code> or <code>GRCh38</code></dd>
</dl>

#### `datasets.csv` data example

| `dataset_name` | `reference_genome` | `file`                                                      |
| :------------- | :----------------- | :---------------------------------------------------------- |
| HG002          | GRCh38             | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002          | GRCh38             | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002          | GRCh38             | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |

### Samples

The `samples.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

#### Data requirements

<dl>
  <dt>sample_name <code>&lt;str&gt;</code></dt>
  <dd>The ID of the sample. this should correspond to the sample ID's provided in the provided <code>.vcf</code> file.
  
  E.g. <code>HG002</code></dd>
  
  <dt>dataset <code>&lt;enum [dataset_name]&gt;</code></dt>
  <dd>The name of the dataset this sample belongs to. This value should correspond to the provided dataset ID listed in <code>datasets.csv</code> 
  
  E.g. <code>1000g</code></dd>
  
  <dt><code>* &lt;str&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis.
  
  E.g. <code>GRCh37</code> or <code>GRCh38</code></dd>
</dl>

#### `samples.csv` data example

| `sample_name` | `dataset` | `SUPER` | `SUB` |
| :------------ | :-------- | :------ | :---- |
| HG002         | HG002     | `EUR`   | `GBR` |
| HG002         | HG003     | `AFR`   | `GWD` |
| HG002         | HG004     | `SAS`   | `GIH` |

### Genomic Locations

The `locations.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

#### Data requirements

<dl>
  <dt>location_name <code>&lt;str&gt;</code></dt>
  <dd>The ID of a gene or, if not a studied gene region, a unique identifier to reference this genomic coordinate window.
  
  E.g. <code>CYP2A6</code></dd>
  
  <dt>chromosome <code>&lt;enum &lt;int [0-24]&gt;&gt;</code></dt>
  <dd>The chromosome number on which the above genomic region can be found.
  
  E.g. <code>19</code></dd>

  <dt>start <code>&lt;int&gt;</code></dt>
  <dd>The start coordinates for the genomic window.
  
  E.g. <code>40842850</code></dd>
  
  <dt>stop <code>&lt;int&gt;</code></dt>
  <dd>The stop coordinates for the genomic window.
  
  E.g. <code>1000g</code></dd>
  
  <dt>strand <code>&lt;enum [-1,1]&gt;</code></dt>
  <dd>The strand on which the genomic region can be found, where <code>1</code> denotes the forward strand and <code>-1</code> denotes the reverse strand.
  
  E.g. <code>-1</code></dd>
</dl>

#### `locations.csv` data example

| `sample_name` | `dataset` | `SUPER`  | `SUB`    |
| :------------ | :-------- | :------- | :------- |
| CYP2A6        | 19        | 40842850 | 40851138 |
| CYP2B6        | 19        | 40988570 | 41021110 |
| UGT2B7        | 4         | 69045214 | 69112987 |

### Transcripts

The `transcripts.csv` file allows you to declare which transcripts you would like to use when performing variant-effect-prediction.

{: .note-title }
> Multiple Transcripts
>
> If more than one transcript is provided for a given genomic region, we will attempt to match the transcripts available in the order that is provided. The first match will be selected, and if no transcripts provided are available, the first available transcript will be selected.

#### Data requirements

<dl>
  <dt>gene_name <code>&lt;enum [str]&gt;</code></dt>
  <dd>The name of the gene a transcript describes. This key should match the gene or region name provided in the <code></code>. 
  
  E.g. <code>HG002</code></dd>
  
  <dt>transcript_id <code>&lt;str&gt;</code></dt>
  <dd>The name of the transcript in question. This value will be used to query the E! Ensembl database when performing variant-effect-prediction. 
  
  E.g. <code>NM_000762.6</code></dd>
</dl>

#### `transcripts.csv` data example

| `gene_name` | `transcript_id`   |
| :---------- | :---------------- |
| CYP2A6      | NM_000762.6       |
| CYP2A6      | ENST00000600495.1 |
| CYP2A6      | ENST00000596719.5 |
| CYP2A6      | ENST00000599960.1 |
| CYP2B6      | NM_000767.5       |
| CYP2B6      | ENST00000593831.1 |
| CYP2B6      | ENST00000598834.2 |
| CYP2B6      | ENST00000597612.1 |
| CYP2B6      | ENST00000594187.1 |
| UGT2B7      | NM_001074.4       |
| UGT2B7      | ENST00000508661.5 |
| UGT2B7      | ENST00000622664.1 |
| UGT2B7      | ENST00000502942.5 |
| UGT2B7      | ENST00000509763.1 |
