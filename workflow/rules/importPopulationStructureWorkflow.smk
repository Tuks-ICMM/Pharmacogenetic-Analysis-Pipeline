
module PopulationStructure:
    snakefile: github("Tuks-ICMM/Population-Structure-Workflow", path="workflow/Snakefile", branch="main")
    config: config

use rule * from PopulationStructure as Population_Structure_*

# TODO: Document this behaviour for technical users
use rule plinkPed from PopulationStructure as Population_Structure_plinkPed with:
    input:
        "results/filter/filterSampleRelatedness.vcf.gz"
# TODO: Document this behaviour for technical users
use rule Plink_PCA from PopulationStructure as Population_Structure_Plink_PCA with:
    input:
        "results/FILTER/ALL_FILTERED.vcf.gz"
# TODO: Document this behaviour for technical users
use rule DAPC from PopulationStructure as Population_Structure_DAPC with:
    input:
        "results/FILTER/ALL_FILTERED.vcf.gz"
