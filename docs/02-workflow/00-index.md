---
title: Workflow
permalink: workflow
layout: page
nav_order: 2
has_children: true
---

# Workflow
{: .no_toc }

A summary of the workflow itself and its analyses, broken down by topic.
{: .fs-6 .fw-300 }

Reference Genome Configuration
{: .label }

---

# Introduction
The <i>{{ site.title }}</i> is a pipeline powered by <a href="https://snakemake.readthedocs.io/" target="_blank">Snakemake</a>, a python-based workflow management package. This project has been created with support for PBS-Torque scheduler environments on Linux servers.

Below is a diagram representing the pipeline flow and steps in the form of a process flow diagram. For reference on the graph syntax (Shape legend), please consult [this guide](https://www.bbc.co.uk/bitesize/guides/znv3rwx/revision/2).

```mermaid
---
title: Pharmacogenetics Analysis Pipeline
---
flowchart TD
    START([Start])
    END([End - Results])

    subgraph ValidateVcfModule
        validateVcf([Validate VCF Workflow])
    end
    subgraph PopulationStructureModule
        PopulationStructureWorkflow([Populaltion Structure Workflow])
    end


    subgraph dataPrep ["Data and Metadata Preparation"]
        %% Use LR to invert axis set by parent to effectively force relative "TB"
        direction LR
        genomeFasta[/"Reference Genome GRCh38 (FASTA)"/]
            subgraph data ["Datasets"]
                datasetFiles1[/"Dataset_1 (.vcf.gz + .tbi)"/]
                datasetFiles2[/"Dataset_2 (.vcf.gz + .tbi)"/]
                datasetFiles3[/"Dataset_... (.vcf.gz + .tbi)"/]
            end
            subgraph metadata ["Analysis metadata"]
                %% Use LR to invert axis set by parent to effectively force relative "TB"
                direction TB

                datasetMeta[/"Datasets metadata (CSV)"/]
                locationMeta[/"Genomic location metadata (CSV)"/]
                sampleMeta[/"Sample metadata (CSV)"/]
                transcriptMeta[/"Transcript metadata (CSV)"/]
            end

    end
    START --> dataPrep
    dataPrep --> validateVcf

    validateVcf --> ifMultipleVcfs

    subgraph prep [Process]
        direction BT

        subgraph multipleVcfProtocol [Multiple dataset protocol]
            mergeVcf[Merge VCFs]
        end


        ifMultipleVcfs{If multiple datasets provided} --> |yes| multipleVcfProtocol
        ifMultipleVcfs{If multiple datasets provided} --> |no| refFromFasta

        refFromFasta[[refFromfasta: Check reference alleles against provided reference genome]]
        
        multipleVcfProtocol --> refFromFasta

        refFromFasta --> chrFilter[[chrFilter: Filter out non-standard chromosomes]]

        chrFilter --> sampleSubset[[sampleSubset: Subset samples to labeled samples in metadata files]]
        
        sampleSubset --> FILTER[[FILTER: Filter variants and samples with 100% missingness & prune overrelated samples]]

        FILTER --> TRIM_AND_NAME[[TRIM_AND_NAME: Trim dataset to genomic regions of interest]]

        TRANSPILE_CLUSTERS[[TRANSPILE_CLUSTERS: Transpile cluster ownership from sample cluster assignment into input format]]

        TRIM_AND_NAME & TRANSPILE_CLUSTERS --> PLINK[[PLINK: Perform frequency analysis]]
 
    end 
    FILTER --> PopulationStructureWorkflow
    PopulationStructureWorkflow --> END
    PLINK --> END
```