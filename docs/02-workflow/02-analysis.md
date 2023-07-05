---
title: Analysis
permalink: workflow/analysis
layout: page
nav_order: 2
parent: Workflow
---

# Configuration
{: .no_toc }

A breakdown of the process used in this workflow and how it has been implemented.
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

```mermaid
---
title: Pharmacogenetics Analysis
---
flowchart TD
subgraph ValidateVcfModule
    validateVcf([Validate VCF Workflow])
end
subgraph PopulationStructureModule
    PopulationStructureWorkflow([Populaltion Structure Workflow])
end

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




## `refFromFasta`
```mermaid
flowchart TD
refFromFasta[[refFromfasta: Check reference alleles against provided reference genome]]
```
This rule is responsible for checking each loci and comparing its listed reference to that provided in the reference genome.

## `chrFilter`
```mermaid
flowchart TD
chrFilter[[chrFilter: Filter out non-standard chromosomes]]
```

## `sampleSubset`
```mermaid
flowchart TD
sampleSubset[[sampleSubset: Subset samples to labeled samples in metadata files]]
```

## `FILTER`
```mermaid
flowchart TD
FILTER[[FILTER: Filter variants and samples with 100% missingness & prune overrelated samples]]
```

## `TRIM_AND_NAME`
```mermaid
flowchart TD
TRIM_AND_NAME[[TRIM_AND_NAME: Trim dataset to genomic regions of interest]]
```

## `TRANSPILE_CLUSTERS`
```mermaid
flowchart TD
TRANSPILE_CLUSTERS[[TRANSPILE_CLUSTERS: Transpile cluster ownership from sample cluster assignment into input format]]
```

## `PLINK`
```mermaid
flowchart TD
PLINK[[PLINK: Perform frequency analysis]]
```