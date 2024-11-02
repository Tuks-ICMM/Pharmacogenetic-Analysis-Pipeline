---
title: Configuration & Data
permalink: workflow/configuration-and-data
layout: page
nav_order: 1
has_children: false
parent: Workflow
---

# Configuration & Data
{: .no_toc}

A summary of the required data and input files needed to perform an analysis.
{: .fs-6 .fw-300 }

<details markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

This page describes the information needed to run the _{{ site.title }}_. Below we guide users through the system used to declare an analysis manifest, and all associated metadata files. For more information, please consult the relevant section below which contains more specific guidance, discussions and technical documentation.

## Overview

This workflow makes use of an analysis manifest to encapsulate all analysis variables used. This manifest file collects and connects the metadata for your samples, datasets, and relevant reference resources (Reference Genomes, etc) together. Doing so allows the workflow to programmatically access clusters through sample annotations, which is required in order to produce cluster-level reports.

<details markdown="block" open>
  <summary>
    Input Data Infographic
  </summary>
  {: .text-delta }

{% raw %}
```mermaid
---
title: Input filemap
config:
    flowchart:
        defaultRenderer: elk
    elk:
        nodePlacementStrategy: BRANDES_KOEPF
---
flowchart TB
  subgraph input [input files]
      subgraph data [Datasets]
          datasetFile1{{<b>Dataset file</b><br><code>GnomAD_Chr1.vcf.gz</code>}}
          datasetFile2{{<b>Dataset file</b><br><code>GnomAD_Chr2.vcf.gz</code>}}
          datasetFileN{{<b>Dataset file</b><br><code>GnomAD_ChrN...vcf.gz</code>}}
      end

      subgraph metadata [Analysis Metadata]
          locationMeta{{<b>Coordinates for study</b><br><code>locations.csv</code>}}
          sampleMeta{{<b>Sample metadata</b><br><code>samples.csv</code>}}
          datasetMeta{{<b>Data files to incude</b><br><code>datasets.csv</code>}}
          transcriptMeta{{<b>Transcript preferences</b><br><code>transcripts.csv</code>}}
      end
  end
  subgraph resources [<code>resources/</code>]
      reference_genome{{Reference Genome <br> <code>resources/genome_version_name.fa</code>}}
  end
  subgraph config [<code>config/</code>]
    configuration{{<b>Analysis configuration</b> <br><code>config/configuration.json</code>}}
  end

  vcf_validation_workflow[\VCF Validation Workflow/]
  click vcf_validation_workflow href "https://tuks-icmm.github.io/VCF-Validation-Workflow/workflow/methodology" _blank

  pharmacogenetic_analysis_workflow[\Pharmacogenetics Analysis Workflow/]
  click pharmacogenetic_analysis_workflow href "/workflow/methodology" _blank

  population_structure_workflow[\Population structure Workflow/]
  click population_structure_workflow href "https://tuks-icmm.github.io/Population-Structure-Workflow/workflow/methodology" _blank

  datasetMeta -...-o|Referenced in| reference_genome

  metadata -.-o|Describes| data

  input --> vcf_validation_workflow
  config ----> vcf_validation_workflow
  resources ----> vcf_validation_workflow

  vcf_validation_workflow --> pharmacogenetic_analysis_workflow

  pharmacogenetic_analysis_workflow --> population_structure_workflow
  pharmacogenetic_analysis_workflow --> results

  population_structure_workflow --> results

  results(((Results)))
```
{% endraw %}

</details>



## Input Data

This workflow is designed to work on variant-call-format files (<code>.vcf</code> file extension). The latest version of the VCF specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.3.pdf).



### Compression and Indexing

This workflow can accept uncompressed VCF files, however this workflow will compress and index the data during handling for performance reasons. If possible, please provide them in compressed and index form.

{: .highlight-title}
> This workflow uses an alternative handling file format
>
> Since a large portion of this workflow makes use of Plink-2, this workflow is configured to convert the files into Plink-2's [binary fileset](https://www.cog-genomics.org/plink/2.0/input#pgen) (<code>.pgen</code>, <code>.pvar</code> and <code>.psam</code> files) offered by [Plink-2](https://www.cog-genomics.org/plink/2.0/), which offers very good processing performance and tracking improvements. We do also perform an additional reference-guided allele verification to remove swapped alleles from data tracked using Plink-1.9 binary files.

## Analysis configuration

To perform an analysis with this workflow, users will need to configure the workflow. This includes providing environment-related information like output locations, as well as analysis settings like reference population selection. This information is all declared and stored using the `config/manifest.json` file.

<h3>The <code>manifest.json</code> file</h3>

This file is responsible for declaring all information relating to the analysis and serves as the central point of contact between the workflow runtime and your input data. It is also used to configure and synchronize any sub-workflows imported internally.

<details markdown="block">
  <summary>
    <code>manifest.json</code> format example
  </summary>
  {: .text-delta }

  <dl>
    <dt><b>input</b> <code>&lt;object&gt;</code></dt>
    <dd>
        <dl>
            <dt><b>datasets</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
            <dd>A list representing the file-path to the dataset metadata file. Should be suitable for use with the python <code>os.path.join()</code> function.</dd>
            <dt><b>locations</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
            <dd>A list representing the file-path to the location metadata file. Should be suitable for use with the python <code>os.path.join()</code> function.</dd>
            <dt><b>samples</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
            <dd>A list representing the file-path to the samples metadata file. Should be suitable for use with the python <code>os.path.join()</code> function.</dd>
            <dt><b>transcripts</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
            <dd>A list representing the file-path to the transcript metadata file. Should be suitable for use with the python <code>os.path.join()</code> function.</dd>
        </dl>
    </dd>
    <dt><b>output</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
    <dd>A list representing a path to a folder where the results of the analysis should be stored. If the folder does not exist, it will be created.</dd>
    <dt><b>resources</b> <code>&lt;Object&gt;</code></dt>
    <dd>
        <dl>
            <dt><b>reference_genomes</b> <code>&lt;Array&lt;Object&gt;&gt;</code></dt>
            <dd>
                This property should contain a list of objects, where each object describes a reference genome available for use, using teh following properties:
                <dl>
                    <dt><b>name</b> <code>&lt;str&gt;</code></dt>
                    <dd>The name of the reference genome. Should correspond to value used in dataset metadata file.</dd>
                    <dt><b>location</b> <code>&lt;Array&lt;str&gt;&gt;</code></dt>
                    <dd>A list representative the file-path to the reference genome. Should be provided in FASTA format.</dd>
                </dl>
            </dd>
        </dl>
    </dd>
    <dt><b>parameters</b> <code>&lt;Object&gt;</code></dt>
    <dd>
        <dl>
            <dt><b>fishers-test</b> <code>&lt;object&gt;</code></dt>
            <dd>
            <dl>
                <dt><i><b>cluster_name</b>*</i> <code>&lt;str&gt;</code></dt>
                <dd>The name of the cluster-level declared in your sample metadata file for which you would like to declare a reference population. This population will be used to conduct pair-wise testing against all remaining populations in the column respectively.</dd>
            </dl>
            </dd>
        </dl>
    </dd>
  </dl>

  ```json
  {
    "input": {
        "datasets": [
            "/",
            "path",
            "to",
            "my",
            "dataset",
            "metadata"
        ],
        "locations": [
            "/",
            "path",
            "to",
            "my",
            "locations",
            "metadata"
        ],
        "samples": [
            "/",
            "path",
            "to",
            "my",
            "samples",
            "metadata"
        ],
        "transcripts": [
            "/",
            "path",
            "to",
            "my",
            "transcripts",
            "metadata"
        ]
    },
    "output": [
        "/",
        "path",
        "to",
        "my",
        "output",
        "location"
    ],
    "resources": {
        "reference_genomes": [
            {
                "name": "",
                "location": [
                    "/",
                    "path",
                    "to",
                    "my",
                    "reference",
                    "genome"
                ]
            }
        ]
    },
    "parameters": {
        "fishers-test": {
            "my_cluster": "my_population_of_interest"
        }
    }
}
  ```
</details>

---

### Metadata

All data and sample metadata is provided in the form of ` .csv` files declared in the `manifest.json` file. These files allow you to declare datasets and provide the necessary information to determine which contig-level files should be used for analysis given the provided genomic coordinates. For convenience, we will assume standard names for the sake of this explanation:

{: .normal }
> This design-pattern of declaring metadata files via the `manifest.json` was chosen specifically to allow users to create and store analysis configurations and metadata alongside data, which often has special storage requirements (e.g. space, access, etc). Through the `manifest.json` file, all other analysis-specific files will be declared and will be accessible. This then only requires that the `manifest.json` file is discoverable under the path `config/manifest.json`, which can be accomplished with a symlink or shortcut, keeping the amount of setup work to a minimum.


#### <code>datasets.csv</code> Metadata

The dataset metadata file allows you to declare information about your datasets to analyze, including the reference genome version and where to locate the files.

{: .highlight }
> Please provide data in the form of multiple <code>*.vcf</code> files split per-contig.

<details markdown="block">
  <summary>
    Format example
  </summary>
  {: .text-delta }

<dl class="def-wide">
  <dt><b>dataset_name</b> <code>&lt;str&gt;</code></dt>
  <dd>The name of the dataset. This value will be used as a universal accessor for that dataset and any information relating to it. This means that any output files will use this value to determine things like filenames, etc. It is also used to connect other metadata to this dataset computationally, E.g. sample-level information.
  
  <br><strong><i>E.g. <code>1000G</code></i></strong></dd>
  
  <dt><b>reference_genome</b> <code>&lt;str&gt;</code></dt>
  <dd>An <code>enum</code> indicating which reference genome version this dataset has been called on.
  
  <br><b><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></b></dd>
  
  <dt><b>file</b> <code>&lt;file_path&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis.
  
  <br><strong><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></strong></dd>
</dl>

| **dataset_name** | **reference_genome** | **file**                                                    |
| :--------------- | :------------------- | :---------------------------------------------------------- |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |
| HG002            | GRCh38               | `/nlustre/users/graeme/PUBLIC/GenomeInABottle/HG002.vcf.gz` |

</details>


#### <code>samples.csv</code> Metadata

The sample metadata file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.

<details markdown="block">
  <summary>
   Format example
  </summary>
  {: .text-delta }

{: .highlight-title }
> Case Sensitive
>
> The following metadata declaration files use _**case-sensitive column names**_.

<dl class="def-wide">
  <dt>sample_name <code>&lt;str&gt;</code></dt>
  <dd>The ID of the sample. this should correspond to the sample ID's provided in the provided <code>.vcf</code> file. 
  
  <br><strong><i>E.g. <code>HG002</code></i></strong></dd>
  
  <dt>dataset <code>&lt;enum [dataset_name]&gt;</code></dt>
  <dd>The name of the dataset this sample belongs to. This value should correspond to the provided dataset ID listed in <code>datasets.csv</code> 
  
  <br><strong><i>E.g. <code>1000g</code></i></strong></dd>
  
  <dt><code>* &lt;str&gt;</code></dt>
  <dd>A file path indicating the location of the dataset to be used in the analysis. Please note that the column names are <b><i><u>case-sensitive</u></i></b>.
  
  <br><strong><i>E.g. <code>GRCh37</code> or <code>GRCh38</code></i></strong></dd>
</dl>

| **sample_name** | **dataset** | **SUPER** | **SUB** |
| :-------------- | :---------- | :-------- | :------ |
| HG002           | HG002       | `EUR`     | `GBR`   |
| HG002           | HG003       | `AFR`     | `GWD`   |
| HG002           | HG004       | `SAS`     | `GIH`   |

</details>


#### <code>locations.csv</code> Metadata

The location metadata file allows you to declare samples and provide the necessary sample-level information for use in this pipeline.


<details markdown="block">
  <summary>
    Format example
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


#### <code>transcripts.csv</code> Metadata

The transcript metadata file allows you to declare which transcripts you would like to use when performing variant-effect-prediction.

During the execution of the _{{ site.title }}_, variant-effect-prediction (VEP) is performed using a publicly accessible VEP query API by E! Ensembl. Currently, the API returns multiple VEP predictions based on any transcripts that are found matching the requested genomic location. Users are able to provide a <code>transcripts.csv</code> input file to declare a list of transcripts per genomic-region they would like to consider for this analysis.

{: .normal-title }
> Transcript IDs
>
>Please use transcripts listed on the [E! Ensembl Database](https://www.ensembl.org/)

{: .normal-title }
> Multiple Transcripts
>
> If more than one transcript is provided for a given genomic region, we will attempt to match the transcripts available in the order that is provided from top to bottom. The first successful VEP transcript match between the users selection and that provided by E! Ensembl will be selected, and if no transcripts provided are available, the first available transcript result will be selected.

<details markdown="block">
  <summary>
   Format example
  </summary>

<dl class="def-wide">
  <dt><b>gene_name</b> <code>&lt;enum [str]&gt;</code></dt>
  <dd>The name of the gene a transcript describes. This key should match the gene or region name provided in the <code>locations.csv</code> file. 
  
  <br><strong><i>E.g. <code>HG002</code></i></strong></dd>
  
  <dt><b>transcript_id</b> <code>&lt;str&gt;</code></dt>
  <dd>The name of the transcript in question. This value will be used to query the E! Ensembl database when performing variant-effect-prediction. 
  
  <br><strong><i>E.g. <code>NM_000762.6</code></i></strong></dd>
</dl>

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

---
