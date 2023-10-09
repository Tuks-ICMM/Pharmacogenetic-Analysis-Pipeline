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
    validateVcf([Validate VCF\nWorkflow])
end
subgraph PopulationStructureModule
    PopulationStructureWorkflow([Populaltion Structure\nWorkflow])
end

validateVcf --> ifMultipleVcfs

subgraph prep [Process]
    direction BT

    subgraph multipleVcfProtocol [Multiple dataset protocol]
        mergeVcf[Merge VCFs]
    end


    ifMultipleVcfs{If multiple\ndatasetsprovided} --> |yes| multipleVcfProtocol
    ifMultipleVcfs --> |no| refFromFasta

    refFromFasta[[refFromfasta:\nCheck reference alleles against\nprovided reference genome]]
    
    multipleVcfProtocol --> refFromFasta

    refFromFasta --> chrFilter[[chrFilter:\nFilter out non-standard\nchromosomes]]

    chrFilter --> sampleSubset[[sampleSubset:\nSubset samples to labeled\nsamples in metadata files]]
    
    sampleSubset --> FILTER[[FILTER:\nFilter variants and samples\nwith 100% missingness & prune\noverrelated samples]]

    FILTER --> TRIM_AND_NAME[[TRIM_AND_NAME:\nTrim dataset to\ngenomic regions of interest]]

    TRANSPILE_CLUSTERS[[TRANSPILE_CLUSTERS:\nTranspile cluster ownership from\nsample cluster assignment into\ninput format]]

    TRIM_AND_NAME & TRANSPILE_CLUSTERS --> PLINK[[PLINK:\nPerform frequency analysis]]

end 
FILTER --> PopulationStructureWorkflow
PopulationStructureWorkflow --> END
PLINK --> END
```


#### Rule Reference

<details markdown="block">
  <summary>
    <code>refFromFasta</code>
  </summary>
  
```mermaid
flowchart TD
refFromFasta[[refFromfasta:\nCheck reference alleles against\nprovided reference genome]]
```
This rule is responsible for checking each loci and comparing its listed reference to that provided in the reference genome.

</details>

<details markdown="block">
  <summary>
    <code>chrFilter</code>
  </summary>

```mermaid
flowchart TD
chrFilter[[chrFilter:\nFilter out non-standard\nchromosomes]]
```

</details>


<details markdown="block">
  <summary>
    <code>sampleSubset</code>
  </summary>

```mermaid
flowchart TD
sampleSubset[[sampleSubset:\nSubset samples to labeled\nsamples in metadata files]]
```

</details>


<details markdown="block">
  <summary>
    <code>FILTER</code>
  </summary>

```mermaid
flowchart TD
FILTER[[FILTER:\nFilter variants and samples\nwith 100% missingness &prune over-related\nsamples]]
```

</details>

<details markdown="block">
  <summary>
    <code>TRIM_AND_NAME</code>
  </summary>
```mermaid
flowchart TD
TRIM_AND_NAME[[TRIM_AND_NAME:\nTrim dataset to genomic regions of interest]]
```
</details>
<details markdown="block">
  <summary>
    <code>TRANSPILE_CLUSTERS</code>
  </summary>
```mermaid
flowchart TD
TRANSPILE_CLUSTERS[[TRANSPILE_CLUSTERS:\nTranspile cluster ownership from\nsample cluster assignment into\ninput format]]
```
</details>
<details markdown="block">
  <summary>
    <code>PLINK</code>
  </summary>
```mermaid
flowchart TD
PLINK[[PLINK:\nPerform frequency analysis]]
```
</details>