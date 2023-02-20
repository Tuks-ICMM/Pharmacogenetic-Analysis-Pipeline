---
layout: page
title: Changelog
permalink: /changelog
nav_order: 5
---

# Changelog
{: .no_toc }

Recent changes to the _{{ site.title }}_.
{: .fs-6 .fw-300 }

Changelog
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
## [v1.0.2-ALPHA](https://github.com/Tuks-ICMM/Vcf-Validation/compare/v1.0.1-ALPHA...v1.0.2-ALPHA) (2023-02-20)

### Features
Changes to `ALL_COLLATE` rule:

- Remove unnecessary variant renaming step. Variant renaming is now done during the `LIFTOVER` rule of the `Vcf-Validation` sub-pipeline only.

Changes to `FILTER` rule:

- Change filter stringency for samples and variants.

- Allow autosomal chromosomes only.

- Remove steps to filter out variants where the reference allele in the `.vcf` file does not match that of the reference genome. This step was unnecessary since variant reference alleles are altered to match that of the reference genome during the `LIFTOVER` rule of the `Vcf-Validation` sub-pipeline.

- Change variant nomenclature specifications for linkage-disequilibrium pruning step.

- Ensure that chromosome output is kept consistent by adding the `--output-chr chr26` flag to all plink commands.


## [v1.0.1-ALPHA](https://github.com/Tuks-ICMM/Vcf-Validation/compare/5a07b1c...v1.0.1-ALPHA) (2023-02-07)

### Features
- Converted 20% missingness tests on variant and sample levels to 100% missingness.
- Included version [v1.1.1-ALPHA](https://github.com/Tuks-ICMM/Vcf-Validation/releases/tag/v1.1.1-ALPHA) of Vcf-Validation external pipeline.


Started tracking active changelog on 24 January 2023.