---
title: Installation
layout: page
permalink: documentation/installation
nav_order: 1
has_children: false
parent: Documentation
---

# Introduction

A brief overview of the workflow, and its proposed use-cases and objectives.
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

## Required Software

### Plink-2

Plink-2 is used throughout this workflow to conduct various data filtering steps, as well as in the calculation of results and generation of reports for downstream processing and graphing. Binaries are available which are ready-to-use upon download on the softwares [homepage](https://www.cog-genomics.org/plink/2.0/).

{: .highlight }
> Due to recent changes relating to Hardy-Wienberg calculations, we are currently using commands that are available trhough the Alpha build.

### Bcftools

This workflow makes use of [Bcftools](https://samtools.github.io/bcftools/bcftools.html) to perform variant normalization and merge operations.

The best way to install this software is to compile a copy for your needs from the source files. Instructions for this can be found [here](https://samtools.github.io/bcftools/howtos/install.html)

### Python

In several places, this workflow makes use of Python scripts to perform data-processing.

## System Requirements

You will need command-line access to execute build and run commands in order to use this workflow.