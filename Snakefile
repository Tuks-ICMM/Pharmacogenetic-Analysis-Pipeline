import os
import pandas as pd

# DEFINE CONFIG FILE FOR SNAKEMAKE:
configfile: "config.json"

# DEFINE CONTEXT-VARIABLES:
finalExtensions=['acount', 'hardy','prune.in', 'prune.out', 'smiss', 'vmiss']
locations=set(config['locations'].keys())
samples=set(config['samples'].keys())
clusters=set(config['cluster']['clusters'])
populations = pd.read_excel("./Clusters.xlsx", engine='openpyxl')
bExtensions=["bed", "bim", "fam"]
tExtensions=["map", "ped"]

# Define a list to store a merge-list of dataset names:
mergeList = list()

# def getReferenceGenome(sample):
#     if config['samples'][wildcards.sample]['refGenome'] == "GRCh37" or config['samples'][wildcards.sample]['refGenome'] == "Hg19":
#         return "/nlustre/users/fourie/H.sapiens/gatk_resource_bundle/2.8/hg19/ucsc.hg19.fasta"
#     else if config['samples'][wildcards.sample]['refGenome'] == "GRCh38":
#         return "/nlustre/users/graeme/pipeline-2020/binaries/hg38.fa.gz"

def directoryExists(path: str):
    if not os.path.exists(path):
                    os.makedirs(path)

onstart:
    directoryExists(".intermediates")
    directoryExists("final")

# BEGIN DEFINING RULES:
rule all:
    """
    Catch-all rule to trigger auto-run of all processes. This process will be fired automatically in absence of explicit process name given by cli-argument.
    """
    input:
        [expand(["final/%s/ALL_{location}.%s.{extension}" % (cluster, population) for population in populations[cluster].unique()], extension=finalExtensions, location=locations)for cluster in ['SUPER', 'SUB']],
        # "final/Admixture/EIGENSOFT.pca",
        # "final/Admixture/EIGENSOFT.plot",
        # "final/Admixture/EIGENSOFT.eval",
        # "final/Admixture/EIGENSOFT.log",
        # "final/Admixture/ADMIXTURE.5.Q",
        # "final/Admixture/ADMIXTURE.5.P"
        ".intermediates/Admixture/ALL.bed",
        ".intermediates/Admixture/ALL.bim",
        ".intermediates/Admixture/ALL.fam",

rule VALIDATE:
    """
    Perform validation of VCF format as well as REF alleles and strip out INFO tags.
    """

    input:
        "input/{sample}.vcf.gz"

    output:
        ".intermediates/PREP/{sample}.vcf.gz",

    params:
        memory="128G"

    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"
    
    run:
        # Remove variant types we cant yet analyse:
        shell("module load gatk-4.0.12.0; gatk SelectVariants  -V {input} --select-type-to-include SNP --select-type-to-include INDEL --select-type-to-exclude MIXED --select-type-to-exclude MNP --select-type-to-exclude SYMBOLIC --exclude-filtered -O .intermediates/PREP/{wildcards.sample}_FILTERED.vcf.gz"),
        # Strip out INFO tags:
        shell("module load bcftools-1.7; bcftools annotate -x INFO -O z -o .intermediates/PREP/{wildcards.sample}_NO_INFO.vcf.gz .intermediates/PREP/{wildcards.sample}_FILTERED.vcf.gz"),
        # Regenerate and verify the VCF header:
        shell("module load picard-2.17.11; java -Xmx{params.memory} -jar $PICARD FixVcfHeader I=.intermediates/PREP/{wildcards.sample}_NO_INFO.vcf.gz O={output}")


rule LIFTOVER:
    """
    Lift Variants onto same Reference build. Otherwise we cant merge them or analyse them in context of each other.
    """
    input:
        ".intermediates/PREP/{sample}.vcf.gz",

    output:
        ".intermediates/LIFTOVER/{sample}.vcf.gz"

    params:
        prefix= lambda wildcards: ".intermediates/LIFTOVER/{sample}_LIFTED".format(sample=wildcards.sample),
        exclusionList= lambda wildcards: ".intermediates/LIFTOVER/{sample}_EXCLUDE.dat".format(sample=wildcards.sample),
        chainFile="binaries/hg19ToHg38.over.chain",
        LiftOver="binaries/liftOverPlink.py",
        rmBadLifts="binaries/rmBadLifts.py",
        ref="binaries/" + config['refGenomes']['GRCh38']

    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"

    run:
        shell("echo 'Determining Liftover requirements now...'")
        if config['samples'][wildcards.sample]['refGenome'] != "GRCh38":
            shell("echo 'Liftover required. All datasets have been mapped to {}'".format(config['samples'][wildcards.sample]['refGenome'])),
            shell("module load liftover"),
            if config['samples'][wildcards.sample]['refGenome'] == "GRCh37" or config['samples'][wildcards.sample]['refGenome'] == "Hg19":
                shell("echo 'Lifting from GRCh37 to GRCh38.'"),
                directoryExists(".intermediates/LIFTOVER")
                shell("module load plink-2; plink2 --vcf .intermediates/PREP/{wildcards.sample}.vcf.gz --set-all-var-ids @:#\$r-\$a --allow-extra-chr --new-id-max-allele-len 40 truncate --chr 1-22 --out .intermediates/LIFTOVER/{wildcards.sample}_PREP --export vcf-4.2 bgz --output-chr chr26"),
                shell("sleep 60; tabix -p vcf .intermediates/LIFTOVER/{wildcards.sample}_PREP.vcf.gz"),
                shell("module load picard-2.17.11; java -Xmx128G -jar $PICARD LiftoverVcf I=.intermediates/LIFTOVER/{wildcards.sample}_PREP.vcf.gz O=.intermediates/LIFTOVER/{wildcards.sample}.vcf.gz C={params.chainFile} REJECT=.intermediates/LIFTOVER/{wildcards.sample}_REJECTED.vcf.gz R={params.ref}"),
        # TODO: Add conditionals for other human reference genome builds
        else:
            print("No liftover required. Dataset {} is already mapped to GRCh38.".format(wildcards.sample)),
            shell("touch .intermediates/LIFTOVER/{wildcards.sample}_EXCLUDE.dat"),
            shell("module load plink-1.9; plink --map input/{wildcards.sample}.map --ped input/{wildcards.sample}.ped --allow-extra-chr --chr 1-22 --recode vcf --keep-allele-order --exclude {params.exclusionList} --out .intermediates/LIFTOVER/{wildcards.sample}"),
        # shell("bgzip .intermediates/LIFTOVER/{wildcards.sample}.vcf"),
        shell("sleep 1m; tabix -f -p vcf .intermediates/LIFTOVER/{wildcards.sample}.vcf.gz"),
        shell("echo '.intermediates/LIFTOVER/{wildcards.sample}.vcf.gz' >> .intermediates/LIFTOVER/merge.list")


rule ALL_COLLATE:
    """
    Collate Datasets together into 1 psudo-dataset for downstream analysis.
    """
    input:
        expand(".intermediates/LIFTOVER/{sample}.vcf.gz", sample=samples)

    output:
        ".intermediates/COLLATE/ALL.vcf.gz",
        ".intermediates/COLLATE/ALL.vcf.gz.tbi"

    params:
        ref="binaries/" + config['refGenomes']['GRCh38']

    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"

    run:
        shell("module load bcftools-1.7; bcftools merge -l .intermediates/LIFTOVER/merge.list -O z -o .intermediates/COLLATE/ALL_PRE.vcf.gz"),
        shell("module load samtools-1.7; tabix -p vcf .intermediates/COLLATE/ALL_PRE.vcf.gz"),
        shell("module load plink-2; plink2 --vcf .intermediates/COLLATE/ALL_PRE.vcf.gz --fa {params.ref} --ref-from-fa force --allow-extra-chr --export vcf bgz --out .intermediates/COLLATE/ALL_REF"),
        shell("module load plink-2; plink2 --vcf .intermediates/COLLATE/ALL_REF.vcf.gz --allow-extra-chr --output-chr chr26 --chr 1-22 --export vcf-4.2 bgz --out .intermediates/COLLATE/ALL"),
        shell("module load samtools-1.7; tabix -p vcf .intermediates/COLLATE/ALL.vcf.gz"),

rule ALL_ANNOTATE:
    """
    Annotate rsID's in psudo-dataset to facilitate down-stream analysis.
    """
    input:
        ".intermediates/COLLATE/ALL.vcf.gz"

    output:
        ".intermediates/ANNOTATE/ALL.vcf.gz"
    
    params:
        refGenome='/apps/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa.gz',
        dbSNP='/nlustre/data/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz'
    
    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"

    run:
        shell("module load bcftools-1.7; bcftools annotate -c ID  -a /nlustre/data/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz -O z -o .intermediates/ANNOTATE/ALL_ANNOTATED.vcf.gz {input}"),
        shell("module load plink-2; plink2 --vcf .intermediates/ANNOTATE/ALL_ANNOTATED.vcf.gz --export vcf-4.2 bgz --out .intermediates/ANNOTATE/ALL"),

rule ADMIXTURE:
    """
    Perform Admixture analysis on the large psudo-dataset (Requires 100 000 minimum variants to distinguish sub-populations and 10 000 to distinguish super-populations.)
    """
    input:
        ".intermediates/ANNOTATE/ALL.vcf.gz"

    output:
        ".intermediates/Admixture/ALL.log",
        ".intermediates/Admixture/ALL.bed",
        ".intermediates/Admixture/ALL.bim",
        ".intermediates/Admixture/ALL.fam",
        "final/Admixture/EIGENSOFT.pca",
        "final/Admixture/EIGENSOFT.plot",
        "final/Admixture/EIGENSOFT.eval",
        "final/Admixture/EIGENSOFT.log",
        "final/Admixture/ADMIXTURE.5.Q",
        "final/Admixture/ADMIXTURE.5.P"

    params:
        path = ".intermediates/Admixture/",
        finalPath = 'final/Admixture/',
        admixtureAssumption = "5"
    
    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"

    run:
        shell("module load bcftools-1.7; bcftools view -O z -o {params.path}FILTERED.vcf.gz -m2 -M2 -v {input}"),
        shell("module load plink-2; plink2 --vcf {params.path}FILTERED.vcf.gz --thin-count 200000 --set-missing-var-ids @_# --make-bed --out {params.path}THINNED"),
        shell("module load admixture-1.3.0; admixture {params.path}THINNED.bed {params.admixtureAssumption}"),
        shell("mkdir {params.finalPath}"),
        shell("cp {params.path}THINNED.{params.admixtureAssumption}.P {params.finalPath}/ADMIXTURE.{params.admixtureAssumption}.P"),
        shell("cp {params.path}THINNED.{params.admixtureAssumption}.Q {params.finalPath}/ADMIXTURE.{params.admixtureAssumption}.Q"),
        shell("mv {params.path}THINNED.bim {params.path}.pedsnp"),
        shell("mv {params.path}THINNED.fam {params.path}.pedind"),
        shell("module load eigensoft; smartpca -i {params.path}THINNED.bed -a {params.path}THINNED.pedsnp -b {params.path}THINNED.pedind -o {params.finalPath}/EIGENSOFT.pca -p {params.finalPath}/EIGENSOFT.plot -e {params.finalPath}/EIGENSOFT.eval -l {params.finalPath}/EIGENSOFT.log")


rule TRIM_AND_NAME:
    """
    Trim the whole-genome psudo-datasets down to several regions of interest for Variant analysis and Variant effect prediction.
    """
    input:
        ".intermediates/ANNOTATE/ALL.vcf.gz"

    output:
        ".intermediates/TRIM/ALL_{location}_TRIMMED.vcf.gz"

    params:
        fromBP = lambda wildcards: config["locations"][wildcards.location]["GRCh38"]["from"],
        toBP = lambda wildcards: config["locations"][wildcards.location]["GRCh38"]["to"],
        chr = lambda wildcards: config["locations"][wildcards.location]["GRCh38"]["chromosome"]
    
    resources:
        cpus=10,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    run:
        shell("module load plink-2; plink2 --vcf {input} --from-bp {params.fromBP} --to-bp {params.toBP} --chr {params.chr} --export vcf-4.2 bgz --out .intermediates/TRIM/ALL_{wildcards.location}_TRIMMED"),


rule ALL_FILTER:
    """
    Filter out individuals missing 100% of their variant information (Safety Check).
    """
    input:
        ".intermediates/TRIM/ALL_{location}_TRIMMED.vcf.gz"

    output:
        "final/ALL_{location}.vcf.gz",
        ".intermediates/FILTER/ALL_{location}_FILTERED.log"
    
    resources:
        cpus=10,
        nodes=1,
        queue="short",
        walltime="00:30:00"

    run:
        shell("module load plink-2; plink2 --vcf {input} --mind 1 --output-chr chr26 --export vcf-4.2 bgz --out .intermediates/FILTER/ALL_{wildcards.location}_FILTERED"),
        shell("module load bcftools-1.7; bcftools norm -m - -o final/ALL_{wildcards.location}.vcf.gz -O z .intermediates/FILTER/ALL_{wildcards.location}_FILTERED.vcf.gz")     

rule TRANSPILE_CLUSTERS:
    """
    Convert Cluster information given in the config file into PLINK-2.0 suitable format.
    """
    output:
        ".intermediates/REFERENCE/cluster_{cluster}.txt"

    resources:
        cpus=15,
        nodes=1,
        queue="normal",
        walltime="30:00:00"
    
    run:
        cluster = pd.read_excel(config['cluster']["file"])
        cluster['FID'] = cluster['ID']
        cluster[['ID', 'FID', wildcards.cluster]].to_csv(".intermediates/REFERENCE/cluster_{}.txt".format(wildcards.cluster), sep='\t', index=False)



rule ALL_ANALYZE:
    """
    Perform Frequency analysis on super populations.
    """
    input:
        expand(".intermediates/REFERENCE/cluster_{cluster}.txt", cluster=clusters),
        vcf="final/ALL_{location}.vcf.gz"
    
    output:
        expand("final/{{cluster}}/ALL_{{location}}.{{subCluster}}.{extension}", extension=finalExtensions,)

    params:
        prefix = 'ALL_{location}'
    
    resources:
        cpus=15,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    run:
        for cluster in clusters:
            shell("module load plink-2; plink2 --vcf {input.vcf} --freq counts --export vcf-4.2 bgz --out final/{cluster}/{params.prefix}"),
            shell("module load plink-2; plink2 --vcf {input.vcf} --double-id --within .intermediates/REFERENCE/cluster_{cluster}.txt populations --freq counts --missing --indep-pairwise 50 5 .05 --hardy midp --loop-cats populations --out final/{cluster}/{params.prefix}"),

# Add in VEP API calls
rule ALL_VEP:
    input:
        vcf=".intermediates/FILTER/ALL_{location}_FILTERED.vcf.gz"

    output:
        excel="final/{location}.xlsx"

    params:
        transcript_id = lambda wildcards: config["locations"][wildcards.location]["GRCh38"]["transcript_id"]

    resources:
        cpus=15,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    script:
        "scripts/VEP.py"



["final/%s/ALL_{{location}}.{subCluster}.{extension}" % cluster for cluster in ['SUPER', 'SUB']]