---
title: Environment
permalink: configuration/environment
layout: page
parent: Workflow configuration
nav_order: 1
has_children: true
---

# Environments
{: .no_toc }

Supported batch schedulers and how to config
{: .fs-6 .fw-300 }

---

Snakemake supports the use of profiles, used to describe the set of commands used to interact with your batch scheduler. A profile is provided as a set of scripts which provide template bach scheduler commands.

{: note-title }
> Warning
>
> Without this profile definition, snakemake will not be able to interact with the batch scheduler directly and will not be able to run rules in parallel.


A profile for the PBS/Torque batch scheduler is available. To create your own `config` file for another scheduler, you can reference [snakemakes documentation]() to learn how. They also provide a link to template scripts for a wider variety of bach schedulers which are sufficient for most uses.
<!-- TODO: Provide link -->
