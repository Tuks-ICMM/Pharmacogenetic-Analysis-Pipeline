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
        validateVcf([Validate VCF\nWorkflow])
    end
    subgraph PopulationStructureModule
        PopulationStructureWorkflow([Populaltion Structure\nWorkflow])
    end


    subgraph dataPrep ["Data and Metadata Preparation"]
        %% Use LR to invert axis set by parent to effectively force relative "TB"
        direction LR
        genomeFasta[/"Reference Genome\nGRCh38 (FASTA)"/]
            subgraph data ["Datasets"]
                datasetFiles1[/"Dataset_1\n(.vcf.gz + .tbi)"/]
                datasetFiles2[/"Dataset_2\n(.vcf.gz + .tbi)"/]
                datasetFiles3[/"Dataset_n...\n(.vcf.gz + .tbi)"/]
            end
            subgraph metadata ["Analysis metadata"]
                %% Use LR to invert axis set by parent to effectively force relative "TB"
                direction TB

                datasetMeta[/"Datasets metadata\n(CSV)"/]
                locationMeta[/"Genomic location\nmetadata (CSV)"/]
                sampleMeta[/"Sample metadata\n(CSV)"/]
                transcriptMeta[/"Transcript metadata\n(CSV)"/]
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


        ifMultipleVcfs{If multiple\ndatasets provided} --> |yes| multipleVcfProtocol
        ifMultipleVcfs --> |no| refFromFasta

        refFromFasta[[refFromfasta:\nCheck reference alleles against\nprovided reference genome]]
        
        multipleVcfProtocol --> refFromFasta

        refFromFasta --> chrFilter[[chrFilter:\nFilter out non-standard\nchromosomes]]

        chrFilter --> sampleSubset[[sampleSubset:\nSubset samples to labeled\nsamples in metadata files]]
        
        sampleSubset --> FILTER[[FILTER:\nFilter variants and samples\nwith 100% missingness & prune\noverrelated samples]]

        FILTER --> TRIM_AND_NAME[[TRIM_AND_NAME:\nTrim dataset to genomic\nregions of interest]]

        TRANSPILE_CLUSTERS[[TRANSPILE_CLUSTERS:\nTranspile cluster ownership\nfrom sample cluster assignment\ninto input format]]

        TRIM_AND_NAME & TRANSPILE_CLUSTERS --> PLINK[[PLINK:\nPerform frequency analysis]]
 
    end 
    FILTER --> PopulationStructureWorkflow
    PopulationStructureWorkflow --> END
    PLINK --> END
```