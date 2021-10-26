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
This pipeline is designed to accept variant call format data in the form of `.vcf` files. Due to some of the bioinformatics software used internally, these files are required to be compressed using BG-zip compression, and provided with an accompanying Tabix index file. Both of these peices of software are provided by [Samtools](https://www.htslib.org/doc/), a standard Bioinformatics software. These two files provide a block-level compressed format of your data, and a block index, allowing the software to decompress portions of your file and access spesific entries without having to decompress the entire file.

>This is also just good practice and **should** be a bioinformatics software standard

> Please be advised, BG-zip compression is **not** the same as gzip compression such as that provided by linuxes gzip command. Though the final output is still block-level compression and is operable by both programs, you will need BG-zip compression in order to create a Tabix index.

### Built with:
This has been made using a python-based [domain spesific language (DSL)](https://www.jetbrains.com/mps/concepts/domain-specific-languages/) called [Snakemake](https://snakemake.readthedocs.io/en/stable/) and coded to run on a PBS/Torque environment using the `qsub` command (this is set by the profile folder). As such, it needs to be run on a server with the appropriate binaries and batch scheduling software. 

#### Software: 
Below is a list of software used by this pipeline:
- [PBS/Torque batch scheduler](https://adaptivecomputing.com/cherry-services/torque-resource-manager/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [PLINK-2.0](https://www.cog-genomics.org/plink2)
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
This pipeline uses the standardised folder structure, where the workflow itself is located under the `workflow` folder.
```
.
├── config                  # All config data (PBS Profile, genes, etc)
├── resources               # Commonly used resources (WARNING: DEPRECIATING SOON)
├── results                 # The output of the pipeline
├── workflow                # The entrypoint to the code of the pipeline
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
3. To generate a runtime report, detailing figures produced and performance-related numbers, use the `--report` snakemake flag _(This requires that you have the [`Jinja2`](https://jinja.palletsprojects.com/en/2.11.x/) python package installed.)_. The HTML file produced is completely self-contained and can be shared as needed. You can view it using any web browser such as firefox or Google Chrome, etc.
___
## Roadmap:
See our [Projects tab](/projects) and [Issues tracker](/issues) for a list of proposed features (and known issues).
___
## Versioning:
We use the [SemVer](http://semver.org/) syntax to manage and maintain version numbers. For the versions available, see the [releases on this repository here](https://github.com/SgtPorkChops/SASDGHUB/releases). 

___
## Acknowledgements:
Many thanks to the following individuals who have been instrumental to the success of this project:
<table>
      <a href="https://github.com/G-kodes">
        <td style="text-align:center;">
          <div>
          <img src=https://avatars0.githubusercontent.com/u/25722914?s=100&v=4" width="100" alt="Graeme Ford" />
          </div>
          <h4>Author</h4>
          <hr/>
          <a href="https://www.linkedin.com/in/graeme-ford/" target="_blank">
            <h4><strong>Graeme Ford</strong></h4>
          </a>
          <a href="https://github.com/orgs/Tuks-ICMM/people/G-kodes" target="_blank">
            <h6><italic>G-Kodes</italic></h6>
          </a>
        </td>
    </a>
      <td style="text-align:center;">
        <div>
        <img src="https://www.up.ac.za/media/shared/489/ZP_Images/michael-pepper-message.zp39643.jpg" width="100" alt="Prof. Michael S. Pepper" />
        </div>
          <h4>Supervisor</h4>
          <hr/>
          <a href="https://www.up.ac.za/institute-for-cellular-and-molecular-medicine/article/2019297/professor-michael-s-pepper" target="_blank">
            <h4><strong>Prof. Michael Pepper</strong></h4>
          </a>
      </td>
      <td style="text-align:center;">
        <div>
        <img src="https://avatars.githubusercontent.com/u/3425899?s=96&v=4" width="100" alt="Prof. Fourie Joubert" />
        </div>
          <h4>Co-Supervisor</h4>
          <hr/>
          <a href="https://www.up.ac.za/the-genomics-research-institute/article/1929131/professor-fourie-joubert" target="_blank">
            <h4><strong>Prof. Fourie Joubert</strong></h4>
          </a>
          <a href="https://github.com/orgs/Tuks-ICMM/people/fouriejoubert" target="_blank">
            <h6><italic>fouriejoubert</italic></h6>
          </a>
      </td>
      <td style="text-align:center;">
        <div>
        <img src="https://avatars.githubusercontent.com/u/87174188?s=96&v=4" width="100" alt="Prof. Fourie Joubert" />
        </div>
          <h4>Tester</h4>
          <hr/>
          <a href="https://www.linkedin.com/in/fatima-barmania-a1201238/" target="_blank">
            <h4><strong>Fatima Barmania-Faruk</strong></h4>
          </a>
          <a href="https://github.com/orgs/Tuks-ICMM/people/Fatimabp" target="_blank">
            <h6><italic>Fatimabp</italic></h6>
          </a>
      </td>
      <td style="text-align:center;">
        <div>
        <img src="https://avatars.githubusercontent.com/u/80751008?s=96&v=4" width="100" alt="PMegan Ryder" />
        </div>
          <h4>Tester</h4>
          <hr/>
          <a href="https://www.linkedin.com/in/megan-ryder-b312b0159/" target="_blank">
            <h4><strong>Megan Ryder</strong></h4>
          </a>
          <a href="https://github.com/orgs/Tuks-ICMM/people/Megs47" target="_blank">
            <h6><italic>Megs47</italic></h6>
          </a>
      </td>
      <td style="text-align:center;">
        <div>
        <img src="https://avatars.githubusercontent.com/u/80751008?s=96&v=4" width="100" alt="PMegan Ryder" />
        </div>
          <h4>Tester</h4>
          <hr/>
          <a href="https://www.linkedin.com/in/megan-ryder-b312b0159/" target="_blank">
            <h4><strong>Sarah Turner</strong></h4>
          </a>
          <a href="https://github.com/orgs/Tuks-ICMM/people/sarahsaraht" target="_blank">
            <h6><italic>sarahsaraht</italic></h6>
          </a>
      </td>
</table>


---