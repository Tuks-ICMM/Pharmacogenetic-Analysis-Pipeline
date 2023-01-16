---
title: Configuration
permalink: overview/configuration
layout: page
nav_order: 2
parent: Overview
---

# Configuration

{: .no-toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

## Environment

The _{{ site.title}}_ currently uses the linux-based PBS-Torque scheduling system. A configuration profile is available under `config/PBS-Torque-Profile` should you wish to expand this profile or otherwise customise it.

{: .normal }

> Contributions and collaborations on addittional platforms profiles are more than welcome.

## Setting global configuration

The _{{ site.title}}_ uses a global configuration located in `config/config.json` to record information that is not analysis-specific. This file contains a top-level JSON `object` to record the configuration options below:

### Reference Genomes

It is possible to set a standard list of available reference genomes under the object-id of `reference_genome` in the form of an `array` of `objects`.

<dl class="def-wide">
  <dt>version <code>&lt;str&gt;</code></dt>
  <dd>The version string to be used to access this reference genome in the pipeline input files.
  
  <br><strong><i>E.g. <code>GRCh38</code></i></strong></dd>
  
  <dt>file_path <code>&lt;array [str]&gt;</code></dt>
  <dd>An array containing the decomposed location of the dataset to be used in the analysis. See the note below for additional information.
  
  <br><strong><i>E.g. <code>["/", "reference", "human", "GRCh38.fa.gz"]</code></i></strong></dd>
</dl>

{: .normal }

> We use the built-in python function `os.path` to generate platform-spesific paths. Should youo wish to provide a path from root, you may do so by setting the first element in the array to the drive reference for your OS. \***\*Linux E.g. ["/", ...]\*\***

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

## Environment options
The _{{ site.title }}_ supports several environmental options which are set at the top-level as follows:

### `environment`
This object contains the configuration for email notifications when a job is done.