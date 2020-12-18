# SnakeMake Pharmacogenetics pipeline
---
## Content Guide:
1. [About the Project](#about-the-project)
    1. [Datasets](#datasets)
    1. [Built with](#built-with)
        1. [Software](#software)
        1. [Binaries](#binaries)
    1. [File Structure](#file-structure)
        1. [Files](#files)
        1. [Folders](#folders)
1. [Usage](#usage)
1. [Roadmap](#roadmap)
1. [Versioning](#versioning)
1. [Authors](#authors)
1. [Acknowledgments](#acknowledgements)

---
## About the project:
This is the development repository for a pipeline created to perform frequency analysis on African genetic datasets.

### Datasets:
The development of this pipeline currently has approval to use datasets from the following sources:
- 1000 Genomes (1000g)
- African Genome Variation Project (AGVP)
- South African Human Genome Project (SAHGP)

### Built with:
This has been made using a python-based [domain spesific language (DSL)](https://www.jetbrains.com/mps/concepts/domain-specific-languages/) called [Snakemake](https://snakemake.readthedocs.io/en/stable/) and coded to run on a PBS/Torque environment using the `qsub` command (this is set by the profile folder). As such, it needs to be run on a server with the appropriate binaries and batch scheduling software. 

#### Software: 
Below is a list of software used by this pipeline:
- PBS/Torque batch scheduler
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [PLINK-1.9](https://www.cog-genomics.org/plink2)
- [VCF-Tools](https://vcftools.github.io/index.html)
- [liftOverPlink](https://github.com/sritchie73/liftOverPlink)(Binaries contained within this repo. _**Update at own risk!**_)
  - [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)(Required dependancy for liftOverPlink)
- [e! Ensembl VEP API](https://www.ensembl.org/info/docs/tools/vep/index.html)
#### Binaries:
Below is a list of binary dependancies used in this pipeline.
- Reference Genomes [_(properly compressed with accompanying index and dictionary files)_](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)
    - GRCh 38
    - Addittional genomes as needed based on input data

### File Structure:
```
.
├── input                   # Analysis-ready data.
├── .intermediates          # Files generated during analysis (hidden).
├── scripts                 # Custom scripts used in analysis.
├── logs                    # SnakeMake Log files.
├── figures                 # Figures generated from analysis.
├── final                   # Final results from analysis.
├── profile                 # Profile setup for PBS/Torque server.
├── binaries                # Reference Material (Ref Genome, etc).
└── README.md
```

This project uses the following naming conventions:

#### Files:
All user generated files should be named using under-score naming conventions. Spaces are replaced with an underscore and co capital letters are used.
> E.g. this_is_a_test_example.txt

All Snakemake generated files are all labeled according to `<sample_name>.<file-extension>` format and stored in a folder named according to the process that produced it.
> E.g. intermediates/liftover/1000g.vcf 

#### Folders:
All user generated folders should use camelCase naming conventions, where the first letter of a multi-word name is lower-case and spaces are removed with the initial letter of the following word capitalised.
> E.g. thisIsATestExample

All snakemake generated folders use the following folder structure:
```
.
└── intermediates
    └── <ruleName>
      └── <file_name>.<extension>
      └── <file_name>.<extension>
      └── <file_name>.<extension> 
```
___
## Usage:
1. use the `cd` command to navigate to the root repository directory containing the `Snakefile`. 
2. To start the pipeline and produce the default list of files, simply call `snakemake` on the command line with appropriate arguments. _(E.g. `--profile` and `--cluster-config` flags)_
3. To generate a runtime report, detailing figures produced and performance-related numbers, use the `--report` snakemake flag. The HTML file produced is completely self-contained and can be shared as needed. You can view it using any web browser such as firefox or Google Chrome, etc.
___
## Roadmap:
See our [Projects tab](/projects) and [Issues tracker](/issues) for a list of proposed features (and known issues).
___
## Versioning:
We use the [SemVer](http://semver.org/) syntax to manage and maintain version numbers. For the versions available, see the [releases on this repository here](https://github.com/SgtPorkChops/SASDGHUB/releases). 
___
## Authors:

<table>
  <a href="https://github.com/G-kodes">
      <tr>
        <td style="text-align:center;">
          <img src=https://avatars0.githubusercontent.com/u/25722914?s=100&v=4" width="100" alt="Graeme Ford" />
        </td>
      </tr>
      <tr>
        <td style="text-align: center;">
          <h4><strong>Graeme Ford</strong></h4>
        </td>
      </tr>
    </a>
</table>

___
## Acknowledgements:
Many thanks to the following individuals who have been instrumental to the success of this project:
<table>
  <a href="https://www.up.ac.za/institute-for-cellular-and-molecular-medicine/article/2019297/professor-michael-s-pepper">
      <td style="text-align:center;">
        <div>
        <img src="https://www.up.ac.za/media/shared/489/ZP_Images/michael-pepper-message.zp39643.jpg" width="100" alt="Prof. Michael S. Pepper" />
        </div>
        <h4>
          <strong>Prof. Michael S. Pepper</strong>
        </h4>
      </td>
    </a>
    <a href="https://www.up.ac.za/institute-for-cellular-and-molecular-medicine/article/2019297/professor-michael-s-pepper">
      <td style="text-align:center;">
        <div>
        <img src="https://www.up.ac.za/media/shared/Legacy/sitefiles/image/proffjoubert.jpg" width="100" alt="Prof. Fourie Joubert" />
        </div>
          <h4><strong>Prof. Fourie Joubert</strong></h4>
      </td>
    </a>
</table>


---