---
title: Introduction
permalink: overview/introduction
nav_order: 1
layout: page
parent: Overview
---

# Introduction
The <i>{{ site.title }}</i> is a pipeline powered by <a href="https://snakemake.readthedocs.io/" target="_blank">Snakemake</a>, a python-based workflow management package. This project has been created with support for PBS-Torque scheduler environments on Linux servers.

Below is a diagram representing the pipeline flow and steps in the form of a process flow diagram. For reference on the graph syntax (Shape legend), please consult [this guide](https://www.bbc.co.uk/bitesize/guides/znv3rwx/revision/2).

```mermaid
---
title: Pharmacogenetics Analysis Pipeline
---
flowchart TB
    START([Start])
    END([End - Results])

    subgraph dataPrep ["Data and Metadata Preparation"]
        %% Use LR to invert axis set by parent to effectively force relative "TB"
        direction LR

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
    end
    START --> dataPrep
    dataPrep --> VALIDATE

    subgraph prep [Pipeline Preparation]


        %% Choice START
        VALIDATE{Validate VCF format} --> |No|VCFERR([Error: Invalid VCF])
        VALIDATE --> |Yes|Qliftover{Reference genome version}
        %% Choice END
        Qliftover --> |GRCh37|LIFTOVER[[Liftover]]
        %% Choice START
        LIFTOVER --> COLLATE
        Qliftover --> |GRCh38|COLLATE[[Collate results into psudo-single dataset]]
        %% Choice END

        COLLATE --> ANNOTATE[[Annotate VCF]]
    end
    ANNOTATE --> processing

    subgraph processing [Data Processing]
        TRIM[[Trim VCF to coordinates of interest]]
    end
    TRIM --> FREQ

    subgraph analysis [Data Analysis]
        ANNOTATE --> |Perform admixture analysis|ADMIXTURE[[Admixture analysis]]

        FREQ[[Perform frequency analysis]]
    end
    FREQ --> END
    ADMIXTURE --> END
```