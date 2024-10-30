---
title: Configuration
permalink: documentation/configuration
layout: page
nav_order: 2
has_children: false
parent: Documentation
---

# Configuration
{: .no_toc }

How to configure the workflow for an analysis
{: .fs-6 .fw-300 }

Software
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

## `config` folder

To perform an analysis with this workflow, a few environmental variables need to be declared. Currently, this includes a reference to the output These are all tracked using files housed in the `config` directory at the project-root. The files in this folder contains all the configuration information specific to the environment you are running this pipeline on.

### `config.json` file

The `config.json` file is used to declare all runtime variables that will be needed to run the workflow and its software.

```json
{
    "fishers-test": {
        "my_cluster": "population_1"
    },
    "reference-genomes": [
        {
            "version": "GRCh38",
            "file_path": [
                "/",
                "path",
                "to",
                "my",
                "reference_genome.fa"
            ]
        }
    ],
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