---
title: PBS/Torque batch scheduler profile
permalink: configuration/environment/pbs-torque
layout: page
nav_order: 1
parent: Environment
grand_parent: Workflow configuration
---

# PBS/Torque batch scheduler profile
{: .no_toc }

Snakemake profile for the submission of scripts trhough the PBS/Torque batch scheduler (`qsub` command).
{: .fs-6 .fw-300 }

PBS/Torque
{: .label }

`qsub`
{: .label }

Schedulers
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

Snakemake provides the ability to declare resource declarations for a wide array of bach scheduler resources such as RAM or CPUs.

The provided batch scheduler profile was created from a template [provided by snakemake](). It has however been altered to facilitate declaring batch scheduler resource usage on a per-rule basis.

{: .note-title }
> Parallelization warning
>
> The ability to set granular resource allocations per-rule was required to achieve parallel rule execution. As such, this pipeline contains rule syntax that is non-standard, to facilitate this paralelization.
>
> Finding ways to standardize this for scalability have been noted and will be investigated when resources allow.
<!-- TODO: Provide link -->