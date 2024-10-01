module VCFValidation:
    snakefile: github("Tuks-ICMM/VCF-Validation-Workflow", path="workflow/Snakefile", branch="Parallelize-Contig-processing")
    config: config

use rule * from VCFValidation as VCF_Validation_*

# TODO: Document this behaviour for technical users
use rule stats from VCFValidation as VCF_Validation_stats with:
    input:
        final = "results/PREP/{dataset_name}/{operation}.vcf.gz",
        final_tbi = "results/PREP/{dataset_name}/{operation}.vcf.gz.tbi",
        initial = lambda wildcards: getattr(rules, str(wildcards.operation)).input.vcf,
        initial_tbi = lambda wildcards: getattr(rules, str(wildcards.operation)).input.tbi

# TODO: Document this behaviour for technical users
ruleorder: VCF_Validation_tabix > tabix


# BEGIN DEFINING RULES:
def vcfValidationWorkflowAdapter(extension: str, wildcards) -> list:
    """
    An adapter to generate the correct input list from `VCF Validation Pipeline`. This is required as liftover is optional, making the output files variable.
    """

    print("TYPE: ", type(wildcards))
    merge_list = list()
    merge_list.extend(expand(outputDir("tmp/sort/{contig}/{dataset_name}_sort{extension}"), dataset_name=dataset_name, extension=extension, contig=wildcards.contig))

    # # [FOR] all unqiue (dataset_name and reference_genome) column combinations present
    # for reference_genome, groupby_subset in datasets.set_index(["reference_genome", "dataset_name"]).groupby(level=0):

    #     # [IF] reference genome version is not "GRCh38" AND there are actual matching records...
    #     if reference_genome != "GRCh38" and groupby_subset is not None:

    #         # [FOR] the column in our MultiIndex that contains the dataset_name's in this subset...
    #         for dataset_name in groupby_subset.index.get_level_values("dataset_name"):
                
    #             # [EACH] add liftover request for the DAG
    #             merge_list.extend(expand(outputDir("tmp/liftover/{contig}/{dataset_name}_liftover{extension}"), dataset_name=dataset_name, extension=extension, contig=datasets["contig"].unique().tolist()))
    #     else:

    #         # [FOR] the column in our MultiIndex that contains the dataset_name's in this subset...
    #         for dataset_name in groupby_subset.index.get_level_values("dataset_name"):
    #                 # [EACH] add liftover request for the DAG
    #                 merge_list.extend(expand(outputDir("tmp/sort/{contig}/{dataset_name}_sort{extension}"), dataset_name=dataset_name, extension=extension, contig=wildcards.contig))
    return merge_list
