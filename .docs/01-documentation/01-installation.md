---
title: Installation
layout: page
permalink: documentation/installation
nav_order: 1
has_children: false
parent: Documentation
---

# Installation

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

This workflow is maintained and distributed through the [ICMM GitHub page](https://github.com/Tuks-ICMM). The source-code has been versioned according to publications and relevant checkpoints, and can be downloaded via a git clone command:

```bash
git clone https://github.com/Tuks-ICMM/VCF-Validation-Workflow.git
```

If you would like to clone a specific version, versions are declared and maintained via the repositories [releases page](https://github.com/Tuks-ICMM/VCF-Validation-Workflow/tags). To access a version, users can clone the repository, and perform a <code>checkout</code> command, providing the version to checkout:

```bash
git clone https://github.com/Tuks-ICMM/VCF-Validation-Workflow.git
git checkout tags/<release_version>
```



## Dependencies

This workflow is built and powered using a python-based framework for workflow management. In addition, we make use of several underlying bioinformatics tools and third-party command-line programs to perform some steps in the analysis (For more information, see the [methodology page](/workflow/methodology)):

- Python
- Plink-2
- BCFtools

### Python

In addition to powering the overall workflow management framework, this workflow makes use of Python scripts to perform data-processing, perform network calls (VEP), and perform data-aggregation where needed. To prevent conflicts with python versions, we recommend the use of conda environments.

{: .normal }
> With increased popularity and adoption of Python, many operating systems now include a stripped down python version for internal use. This does pose some logistical issues when a user may want to install a custom, complete installation for direct use. It is possible in many cases to mix and overwrite python versions unless a dedicated strategy or version-management system is used. This can be especially dangerous as depending on the damage to the OS python installation, fixing this issue may require re-installing the OS.

### Plink-2

Plink-2 is a free, open-source genome association analysis cli program which provides functions to perform several different types of bioinformatics analysis at scale in a computationally efficient setting. This program used throughout this workflow to conduct various data filtering steps, as well as to calculate and compile reports used for downstream processing and graphing.

Installation methods include pre-compiled binaries, available and ready-to-use upon download on the softwares [homepage](https://www.cog-genomics.org/plink/2.0/). Once downloaded, the contained executable should be made accessible via the systems environmental `PATH` variable.

{: .normal }
> The [source code](https://github.com/chrchang/plink-ng/tree/master/2.0/) is available, should users need to compile a copy fo the software for their system. The steps to compile this software are beyond the scope of this guide.

{: .highlight }
> Due to recent changes relating to Hardy-Weinberg calculations, we are currently using commands that are available through the Alpha build.

You may test the installation by opening your command-line and executing the following command:

```
plink2 --help
```

### BCFtools

[BCFtools](https://samtools.github.io/bcftools/bcftools.html) is a command-line bioinformatics program which provides a collection of functions to allow for the manipulation of variant call records found in VCF files. It has been used in several places to perform VCF file merging, normalization and other quality-control operations.

The best way to install this software is to compile a copy for your needs from the source files. Instructions for this can be found [here](https://samtools.github.io/bcftools/howtos/install.html)

{: .normal }
> Pre-compiled copies are available for use through other platforms such as package-managers (E.g. Ubuntu's `apt`), however these are typically not the most up-to-date version of the tool. For this reason, we recommend making use of a copy compiled from source for your needs, as this will grantee the latest available version.
