module VCFValidation:
    snakefile: github("Tuks-ICMM/VCF-Validation-Workflow", path="workflow/Snakefile", branch="main")
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