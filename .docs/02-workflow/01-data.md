---
title: Data
permalink: workflow/data
layout: page
nav_order: 1
has_children: false
parent: Workflow
---

# Data
{: .no_toc}

A summary of the required data and input files needed to perform an analysis.
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

This page lists the information needed to run the _{{ site.title }}_. Below we guide users through the system used to declare an analysis manifest, and all associated metadata files. For more information, please consult the relevant section below which contains more specific guidance, discussions and technical documentation.

## Overview

This workflow makes use of an analysis manifest to encapsulate all analysis variables used. This manifest file collects and connects the metadata for your samples, datasets, and relevant reference resources (Reference Genomes, etc) together. Doing so allows the workflow to programmatically access clusters through sample annotations, which is required in order to produce cluster-level reports.

<details markdown="block" open>
  <summary>
    Input Data Infographic
  </summary>
  {: .text-delta }

```mermaid
---
title: Data Requirements
---
{% raw %}
flowchart TB
  subgraph Standard [Resources folder]
      reference_genome{{Reference Genome <br> <i>genome_version_name.fa</i>}}
  end
  subgraph projectSpecific [Input folder]
      subgraph data [Analysis datasets]
          direction TB
          datasetFile1{{<b>Dataset file</b><br><code>input/GnomAD_Chr1.vcf.gz</code>}}
          datasetFile2{{<b>Dataset file</b><br><code>input/GnomAD_Chr2.vcf.gz</code>}}
          datasetFileN{{<b>Dataset file</b><br><code>input/GnomAD_ChrN...vcf.gz</code>}}
      end
      subgraph metadata [Analysis Metadata]
          direction LR

          manifest{{<b>Analysis Manifest</b> <br><code>input/manifest.json</code>}}

          datasetMeta{{<b>Data files to incude</b><br><code>input/datasets.csv</code>}}
          locationMeta{{<b>Coordinates for study</b><br><code>input/locations.csv</code>}}
          sampleMeta{{<b>Sample metadata</b><br><code>input/samples.csv</code>}}
          transcriptMeta{{<b>Transcript preferences</b><br><code>input/transcripts.csv</code>}}
      end
  end

  workflow[\Pharmacogenetics Analysis Workflow/]

  click workflow href "/workflow/methodology" _blank
  
  transcriptMeta & sampleMeta & datasetMeta & locationMeta --> manifest

  data -.-|Referenced in| datasetMeta
  data -.-|Referenced in| sampleMeta
  data -.-|Referenced in| locationMeta

  reference_genome -.-|Referenced in| datasetMeta

  manifest --> workflow

{% endraw %}
```

</details>


## Analysis Datasets

This workflow is designed to work on <i>variant-call-format</i> files (<code>.vcf</code> file extension) files. The latest version of the VCF specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.3.pdf).

### Dataset Subdivisions

The data used in this workflow should be sub-divided into contigs. This reduces unnecessary processing times associated with genomic content that is not relevant to the analysis.

### Dataset Compression and Indexing

Datasets are often quite large in uncompressed form. Users are welcome to compress their data files for additional performance gains. The software used in this workflow supports BGZip-compression.

If you wish to compress your VCF files, please provide the following files as input:

- BGZIP-compressed VCF file (<code>.vcf.gz</code> or <code>vcf.bgz</code>)
- Tabix Index (<code>.vcf.gz.tbi</code> or <code>.vcf.bgz.tbi</code>)

{: .normal }
> <b>Block Compression (BGZIP)</b> is a non-standard type of compression which is not the default compression type used on Windows or MacOS. It is used to compress files in a series of blocks or chunks.
>
> Block-compression alone is simply an alternative compression method to make your data file smaller. In computational biology applications, block-compression is combined with a <b>Tabix Index</b> to record the coordinate coverage/bounds in each compressed block. This allows targeted decompression of spesific regions for analysis, as opposed to having to parse the entire file until the requested coordinates are found.
>
> Both block-compression and tabix indexing are provided as part of [SamTools](http://www.htslib.org/doc/bgzip.html).

## Analysis Metadata

<h3><code>manifest.json</code></h3>

To declare all analysis variables, a <code>manifest.json</code> file should be provided in the <code>input</code> folder. This file is responsible for declaring all information relating to the analysis and serves as the central point of contact between the workflow runtime and your input data.

All metadata is provided in the form of appropriately named ` .csv` files located in the input directory.

{: .normal-title }
> Case sensitivity
>
> The following metadata declaration files use _**case-sensitive column names**_.

<details markdown="block">
  <summary>
    <h4><code>manifest.json</code> format example</h4>
  </summary>

  <dl>
    <dt>fishers-test <code>&lt;object&gt;</code></dt>
    <dd>
      <dl>
        <dt><i>cluster_name</i> <code>&lt;str&gt;</code></dt>
        <dd>The name of the cluster-level declared in your <code>samples.csv</code> annotations for which you would like to declare a reference population for pair-wise testing.</dd>
      </dl>
    </dd>
    <dt>output-dir <code>&lt;Array&lt;Str&gt;&gt;</code></dt>
    <dd>A list representing the file-path for the location at which the workflow should save its output. If the folder does not exist, the workflow will automatically create it.</dd>
  </dl>

  ```json
  {
      "fishers-test": {
          "my_cluster": "my_population_of_interest"
      },
      "output-dir": [
          "/",
          "path",
          "to",
          "my",
          "output",
          "location"
      ]
  }
  ```
</details>

---
<h3><code>datasets.csv</code></h3>

The `datasets.csv` file allows you to declare datasets and provide the necessary information to determine which contig-level files should be used for analysis given the provided genomic coordinates.


<details markdown="block">
  <summary>
    <h4><code>datasets.csv</code> format example</h4>
  </summary>
  {: .text-delta }

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

| **dataset_name** | **reference_genome** | **file**                                                    |
| :--------------- | :------------------- | :---------------------------------------------------------- |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |

</details>

---
<h3><code>samples.csv</code></h3>

The `samples.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

<details markdown="block">
  <summary>
    <h4><code>samples.csv</code> format example</h4>
  </summary>
  {: .text-delta }

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
| **sample_name** | **dataset** | **SUPER** | **SUB** |
| :-------------- | :---------- | :-------- | :------ |
| HG002           | HG002       | `EUR`     | `GBR`   |
| HG002           | HG003       | `AFR`     | `GWD`   |
| HG002           | HG004       | `SAS`     | `GIH`   |

</details>

---
<h3><code>locations.csv</code></h3>

The `locations.csv` file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.


<details markdown="block">
  <summary>
    <h4><code>locations.csv</code> format example</h4>
  </summary>
  {: .text-delta }


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
| **location_name** | **chromosome** | **start** | **stop**  | **strand** |
| :---------------- | :------------- | :-------- | :-------- | :--------- |
| CYP2A6            | 19             | 40842850  | 40851138  | -1         |
| CYP2B6            | 19             | 40988570  | 41021110  | 1          |
| UGT2B7            | 4              | 69045214  | 69112987  | 1          |

</details>

---
### <code>transcripts.csv</code>

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

<details markdown="block">
  <summary>
    <code>transcripts.csv</code> data example
  </summary>

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

</details>



## Reference Genome

Reference Genomes are considered a 'resource', given that they are versioned and in many cases, may be analysis-specific and non-standard. As such, they are housed under the `resources` folder, and need to be provisioned and maintained separately.
