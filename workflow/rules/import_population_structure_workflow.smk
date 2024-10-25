module population_structure:
    snakefile: github("Tuks-ICMM/Population-Structure-Workflow", path="workflow/Snakefile", branch="main")
    config: config

use rule * from population_structure as population_structure_*

# TODO: Document this behaviour for technical users
use rule remove_rare_variants from population_structure as population_structure_remove_rare_variants with:
    input:
        pgen=outputDir("tmp/{contig}/removed_related_samples.pgen"),
        pvar=outputDir("tmp/{contig}/removed_related_samples.pvar.zst"),
        psam=outputDir("tmp/{contig}/removed_related_samples.psam"),
        sample_metadata=outputDir("tmp/formatted_sample_metadata/samples.tsv")

# # TODO: Document this behaviour for technical users
# use rule plinkPed from PopulationStructure as Population_Structure_plinkPed with:
#     input:
#         "results/filter/filterSampleRelatedness.vcf.gz"
# # TODO: Document this behaviour for technical users
# use rule Plink_PCA from PopulationStructure as Population_Structure_Plink_PCA with:
#     input:
#         "results/filter/filterSampleRelatedness.vcf.gz"
# # TODO: Document this behaviour for technical users
# use rule DAPC from PopulationStructure as Population_Structure_DAPC with:
#     input:
#         "results/filter/filterSampleRelatedness.vcf.gz"
