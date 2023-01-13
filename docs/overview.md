---
title: Overview
permalink: /overview
nav_order: 3
---

# Overview

The <i>{{ site.title }}</i> is a pipeline powered by <a href="https://snakemake.readthedocs.io/" target="_blank">Snakemake</a>, a python-based workflow management package. This project has been created with support for PBS-Torque scheduler environments on Linux servers.

Below is a diagram representing the pipeline flow and steps:

```mermaid
flowchart TB

    subgraph prep [Data Preparation]
        ALL --> VALIDATE
        VALIDATE -.->|Different reference Genomes| LIFTOVER
        LIFTOVER --> COLLATE
        VALIDATE -.->|Same Reference Genomes| COLLATE
        COLLATE --> ANNOTATE
    end
    subgraph processing [Data Processing]
        prep --> TRIM
        TRIM --> FILTER
    end

    subgraph admix [Admixture]
    end

    subgraph analysis [Data Analysis]
        processing ---> PLINK
        processing ---> admix
    end

    prep --> processing
```
