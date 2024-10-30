---
title: Usage
permalink: documentation/usage
layout: page
nav_order: 3
has_children: false
parent: Documentation
---

# Usage
{: .no_toc }

A brief guide to using Snakemake to execute an analysis
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

This workflow is powered using [Snakemake](https://snakemake.readthedocs.io/en/stable/), a python-based workflow-management framework. As such, running this workflow will require an understanding of Snakemake and their underlying philosophy around workflow/pipelines. Below is a brief summary of the important points needed to begin running this workflow, however an in-depth review of the Snakemake framework is reccomended.

## DAG

Each <code>rule</code> or <code>checkpoint</code> in this or any connected workflow(s) is defined as a collection of properties, declaring various properties associated with the <code>rule</code>/<code>checkpoint</code>. 

By understanding which <code>rule</code>/<code>checkpoint</code> produces a given output, and which <code>rule</code>/<code>checkpoint</code> consumes that output as its input, Snakemake is able to chain together the <code>rule</code>/<code>checkpoint</code> and determine the order in which these pieces of code should be executed. Snakemake constructs a  directed-acyclic-graph (DAG) at runtime, using the input and output dependencies declared by each <code>rule</code>/<code>checkpoint</code> available in this or any connected workflows to align them and construct a dependency chain.

As a result of this approach, it is theoretically possible to declare complex DAGs where entire segments of the DAG are never run, or are only run under spesific circumstances. It is then not out of the question to end up with a DAG which contains unconnected clouds of rules, which raises several questions regarding teh handling of large, complex trees of rules.

As a result of this, <b>Snakemake works backwards from the requested files when determining which rules need to be run</b>. What this means is that when a Snakemake workflow is executed, requested output files need to be provided to give Snakemake an anchor-point from which to construct the DAG. It is possible to declare default outputs in a rule from the workflow (titled <code>all</code> by and annotated with a flag marking it as default). This <code>all</code> rule provides an overridable list of standard files to generate when the execution command is run without explicit file requests.