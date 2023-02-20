---
title: Configuration
permalink: overview/configuration
layout: page
nav_order: 2
parent: Overview
---

# Configuration
{: .no_toc }

How to set up infrastructure-related settings and provide standard analysis-agnostic data files.
{: .fs-6 .fw-300 }

Reference Genome Configuration
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

## Environment
The _{{ site.title }}_ currently uses the linux-based PBS-Torque scheduling system. A configuration profile is available under `config/PBS-Torque-Profile` should you wish to expand this profile or otherwise customize it.

{: .normal }
> Contributions and collaborations on additional platforms profiles are more than welcome.

## Setting global configuration

The _{{ site.title }}_ uses a global configuration located in `config/config.json` to record information that is not analysis-specific. This file contains a top-level JSON `object` to record the configuration options below:
### Reference Genomes

It is possible to set a standard list of available reference genomes under the object-id of `reference_genome` in the form of an `array` of `objects`.

<dl class="def-wide">
  <dt>version
    <code>&lt;str&gt;</code>
  </dt>
  <dd>The version string to be used to access this reference genome in the pipeline input files.
    <br><strong>
      <i>E.g.
        <code>GRCh38</code>
      </i>
    </strong>
  </dd>

  <dt>file_path
    <code>&lt;array [str]&gt;</code>
  </dt>
  <dd>An array containing the decomposed location of the dataset to be used in the analysis. See the note below for additional information.
    <br><strong>
      <i>E.g.
        <code>["/", "reference", "human", "GRCh38.fa.gz"]</code>
      </i>
    </strong>
  </dd>
</dl>

{: .normal }
> We use the built-in python function `os.path` to generate platform-specific paths. Should you wish to provide a path from root, you may do so by setting the first element in the array to the drive reference for your OS. \***\*Linux E.g. ["/", ...]\*\***

Users may also opt to use compression and indexing for performance gain when using reference genomes. In such cases, only the declaration for the compressed `FASTA` file is required in `config.json`. Block compression (BGZIP) such as that provided by [SamTools](http://www.htslib.org/doc/bgzip.html) can be used to compress these `FASTA` files. In such cases, additional accompanying index and dictionary files will be required to facilitate performant and targeted decompression. These accompanying files will need to be stored alongside this file. 

These files include:
-  a `.fa.gz.gzi` file (Can be generated during compression) ([Samtools](http://www.htslib.org/doc/bgzip.html))
-  a `.dict` file ([Samtools](http://www.htslib.org/doc/samtools-dict.html))
-  a `.fa.gz.fai` file ([Samtools](http://www.htslib.org/doc/samtools-faidx.html))

**Example:**

```json
{
  "reference_genome": [
    {
      "version": "GRCh38",
      "file_path": ["/", "reference", "human", "GRCh38.fa.gz"]
    },
    {
      "version": "GRCh37",
      "file_path": ["/", "reference", "human", "GRCh37.fa.gz"]
    }
  ]
}
```

---
## Environment options

The _{{ site.title }}_ supports several environmental options which are set at the top-level as follows:

### `environment`

This object contains the configuration for all infrastructure-related configurations.

---
#### `email`
PBS-Torque queue scheduling systems currently supports notification emails. To use this feature, you may provide the following information:

<dl class="def-wide">
  <dt>email <code>&lt;str [Email]&gt;</code></dt>
  <dd>An email address to which the notification should be sent.</dd>

  <dt>conditions <code>[ &lt;enum ['o', 'e']&gt; ]</code></dt>
  <dd>An array of mail-options which indicates when you should receive a notification email for this pipeline execution. <code>a</code> indicates mail should be sent when job is aborted, <code>b</code> indicates mail should be sent when job begins and <code>e</code> indicates mail should be sent when job terminates.</dd>
</dl>

##### `email` example

---
#### `working-directory`
This `object` property is used to denote the current working directory for internal reference purposes.

##### `working-directory` example
```json
{
  "working-directory": "/my/path/"
}
```

#### `queues`
The PBS-Torque batch scheduling system supports the creation of generic resources available for request by users. The _{{ site.title }}_ has been designed to take advantage of this queue system and split each step in the workflow into a separate job submission. This allows us to parallelize as much of the analysis as possible, providing performance bonuses. Custom installations, however, may contain custom queue definitions which are managed on a per-installation manner and typically, are managed by sysAdmins.

You may use the `queue` key to provide an array of objects declaring each type of PBS-Torque queue and its resource availabilities. The properties in this `queue` object correspond to PBS-Torque resource restrictions which should be provided by the administrator/s of your cluster.

{: .note-title }
> Custom core and node selections
>
> In some cases, users might want to run some jobs on multiple nodes and some on single nodes. To support this, you may declare the same underlying queue multiple times with a different `queue` key in teh config file and create multiple versions of the same underlying hardware queue.

{: .note }
> It is recommended that you submit the `all` rule with the longest 

<dl>
  <dt>queue <code>&lt;str&gt;</code></dt>
  <dd>The name of the queue.</dd>
  
  <dt>walltime <code>&lt;str&gt;</code></dt>
  <dd>The maximum walltime jobs on this queue are permitted to execute in a HH:MM:SS format.
  <br><strong><i>E.g. "900:00:00" = 37.5 days</i></strong></dd>

  <dt>memory <code>&lt;str&gt;</code></dt>
  <dd>The amount of RAM available on this queue.
  
  <br><strong><i>E.g. 128G</i></strong></dd>

  <dt>cores <code>&lt;str&gt;</code></dt>
  <dd>The number of cores available on this queue.
  
  <br><strong><i>E.g. 10</i></strong></dd>

  <dt>nodes <code>&lt;str&gt;</code></dt>
  <dd>The number of nodes available in this queue.
  
  <br><strong><i>E.g. 1</i></strong></dd>

  <dt>rules <code>&lt;array [&lt;str&gt;]&gt;</code></dt>
  <dd>An array of rules this rule should be used for. For a reference of rules, please reference the rules list included in teh example below.</dd>
</dl>


##### `queues` example
```json
{
  "queues": [
    {
      "queue": "long",
      "walltime": "900:00:00",
      "memory": "128G",
      "cores": "10",
      "nodes": "1",
      "rules": [
        "all",
        "VALIDATE",
        "LIFTOVER",
        "COLLATE",
        "ALL_COLLATE",
        "ANNOTATE",
        "ADMIXTURE",
        "TRIM_AND_NAME",
        "FILTER",
        "TRANSPILE_CLUSTERS",
        "PLINK"
      ]
    }
  ]
}
```