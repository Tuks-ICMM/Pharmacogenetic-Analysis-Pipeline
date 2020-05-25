import os

configfile: "config.json"

def get_sec(time_str):
    """Get Seconds from time."""
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

finalExtensions=['frq.strat', 'fst', 'imiss', 'ld', 'missing.hap', 'hwe', 'het', 'ibc']
locations=list(config['locations'].keys())
samples=list(config['samples'].keys())
superPop=set(config['clusters']['SUPER'])
subPop=set(config['clusters']['SUB'])
bExtensions=["bed", "bim", "fam"]
tExtensions=["map", "ped"]

rule all:
    input:
        expand("Final/SUB/ALL_{location}_SUB.{extension}", extension=finalExtensions, location=locations),
        expand("Final/SUPER/ALL_{location}_SUPER.{extension}", extension=finalExtensions, location=locations)


rule TRIM_AND_NAME:
    input:
        expand("rawData/{{sample}}.{extension}", extension=bExtensions)

    output:
        expand("Intermediates/TRIM/{{sample}}_{{location}}_READY.{extension}", extension=tExtensions)

    params:
        fromBP = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["from"],
        toBP = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["to"],
        chr = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["chromosome"]

    shell:
        """
        module load plink-1.9
        module load plink2
        plink --bfile rawData/{wildcards.sample} --chr {params.chr} --set-missing-var-ids @_# --make-bed --keep-allele-order --from-bp {params.fromBP} --to-bp {params.toBP} --out Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_TRIMMED
        plink2 --bfile Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_TRIMMED --set-all-var-ids @:#\$r-\$a --new-id-max-allele-len 40 truncate  --make-bed --out Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_NAMED
        plink --bfile Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_NAMED --keep-allele-order --recode --out Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_READY
        """


#/********* LiftOverPlink (LIFTOVER) *********/

rule LIFTOVER:
    input:
        expand("Intermediates/TRIM/{{sample}}_{{location}}_READY.{extension}", extension=tExtensions)

    output:
        outMap="Intermediates/LIFTOVER/{sample}_{location}_LIFTED.map",
        outPed="Intermediates/LIFTOVER/{sample}_{location}_LIFTED.ped",
        exclusion="Intermediates/LIFTOVER/{sample}_{location}_EXCLUDE.dat"

    params:
        prefix=lambda wildcards: "Intermediates/LIFTOVER/{}_{}_PRE_FILTER".format(wildcards.sample, wildcards.location),
        prefix2=lambda wildcards: "Intermediates/LIFTOVER/{}_{}_POST_FILTER".format(wildcards.sample, wildcards.location),
        prefix3=lambda wildcards: "Intermediates/LIFTOVER/{}_{}_EXCLUDE.dat".format(wildcards.sample, wildcards.location),
        chainFile="Binaries/hg19ToHg38.over.chain",
        LiftOver="Binaries/liftOverPlink.py",
        rmBadLifts="Binaries/rmBadLifts.py",
        sexes="rawData/1000g.sexes"

    run:
        shell("echo 'Determining Liftover requirements now...'")
        if config['samples'][wildcards.sample]['refGenome'] != "GRCh38":
            shell("echo 'Liftover required. Dataset {} is mapped to {}'".format(wildcards.sample, config['samples'][wildcards.sample]['refGenome'])),
            shell("module load liftover"),
            if config['samples'][wildcards.sample]['refGenome'] == "GRCh37" or config['samples'][wildcards.sample]['refGenome'] == "Hg19":
                shell("echo 'Lifting from GRCh37 to GRCh38.'"),
                shell("python {params.LiftOver} -e /apps/liftover/liftOver -m Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_READY.map -p Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_READY.ped -o {params.prefix} -c {params.chainFile}"),
                shell("echo 'liftOver complete. Removing Bad lifts.'"),
                shell("python {params.rmBadLifts} --map {params.prefix}.map --out {output.outMap} --log {params.prefix2}.log"),
                shell("echo 'Bad lifts removed.'"),
                shell("cut -f 2 {params.prefix2}.log > {params.prefix3}"),
                shell("chmod a+rwx {params.prefix3}"),
                shell("cut -f 4 {params.prefix}.unlifted | sed '/^#/d' >> {params.prefix3}"),
                shell("less {params.prefix3} | {output.exclusion}"),
                shell("echo 'Exclusion list generated'")
                if config['samples'][wildcards.sample]['requiresSexing']:
                    shell("echo 'config.json file indicates sexing is required for {} dataset.'".format(wildcards.sample)),
                    shell("awk 'NR==FNR{{a[$1]=$4;next}}{{ if (a[$1]!= 'NULL') {{ if (a[$1] == 'female') {{ $5=2; print }} else if (a[$1] == 'male') {{ $5=1; print }} else {{$5=0; print}}}}}}' {params.sexes} {params.prefix}.ped > Intermediates/LIFTOVER/{wildcards.sample}_{wildcards.location}_LIFTED.ped"),
                    shell("echo 'Sexing complete'")
                else:
                    shell("less {params.prefix}.ped > Intermediates/LIFTOVER/{wildcards.sample}_{wildcards.location}_LIFTED.ped")
        # ToDo: Add conditionals for other human reference genome builds
        else:
            print("No liftover required. Dataset {} is already mapped to GRCh38.".format(wildcards.sample)),
            shell("cp Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_READY.map Intermediates/LIFTOVER/{wildcards.sample}_{wildcards.location}_LIFTED.map"),
            shell("cp Intermediates/TRIM/{wildcards.sample}_{wildcards.location}_READY.ped Intermediates/LIFTOVER/{wildcards.sample}_{wildcards.location}_LIFTED.ped"),
            shell("touch Intermediates/LIFTOVER/{wildcards.sample}_{wildcards.location}_EXCLUDE.dat")


rule CLEAN:
    input:
        mapFile="Intermediates/LIFTOVER/{sample}_{location}_LIFTED.map",
        pedFile="Intermediates/LIFTOVER/{sample}_{location}_LIFTED.ped",
        exclude="Intermediates/LIFTOVER/{sample}_{location}_EXCLUDE.dat"

    output:
        expand("Intermediates/CLEAN/{{sample}}_{{location}}_CLEANED.{extension}", extension=['bed', 'bim', 'fam'])

    shell:
        """
        module load plink-1.9
        plink --map {input.mapFile} --ped {input.pedFile} --out Intermediates/CLEAN/{wildcards.sample}_{wildcards.location}_CLEANED --make-bed --keep-allele-order --exclude {input.exclude}
        """

# rule FST:
#     label 'CYP2A6Analysis'
#     label 'normalQueue'

#     input:
#     set file(BEDFile), file(BIMFile), file(FAMFile) from G_1_cleaned_1
#     file SubPop from G_1_PLINK_SUBPOP
    
#     output:
    
#     script:
#     """
#     module load plink-1.9
#     plink --bed {BEDFile} --bim {BIMFile} --fam {FAMFile} --within {SubPop} --fst --out 5-G_1_SUBFST
#     """
def get_prev(input):
    print(input.filename)
    new = ""
    for i in input.filename:
        if i == "-":
            new += "-"
        else:
            pass

    large = max(new.split(), key=len)
    print(input.filename.split(large))
    extensions = [".bed", ".bim", ".fam"]
    return ("Intermediates/" + x + y for x in input.filename.split(large) for y in extensions)

def getFinalName(datasets):
    hold = ""
    batch_a = list()
    batch_b = list()
    final = list(datasets)


    c = 0
    while len(final) > 1:
        if len(final) % 2 > 0:
            hold = final[-1]
            del final[-1]
            batch_a = final[::2]
            batch_b = final[1::2]
            final = [hold]
        else:
            batch_a = final[::2]
            batch_b = final[1::2]
            final = []
        
        c += 1
        for i,j in zip(batch_a, batch_b):
            # Merge files
            final.append(f"{i}{'-' * c}{j}")
    return str(final[0])

rule ALL_COLLATE:
    input:
        #get_prev
        expand("Intermediates/CLEAN/{sample}_{{location}}_CLEANED.{extension}", extension=bExtensions, sample=samples),
        #supPopClusters="rawData/superPopCluster",
        #subPopClusters="rawData/subPopCluster",
    
    output:
        expand("Intermediates/COLLATE_{{location}}/ALL_{{location}}.{extension}", extension=bExtensions)
        #filter(lambda fn: all(e in fn for e in cond2), fns)
        #expand("Intermediates/ALL_{{location}}.{extension}", extension=['bed', 'bim', 'fam'])

    params:
        prefix = "ALL_1_PRE_COLLATE",
        prefix2 = "ALL_1_MERGE",
        prefix3 = "ALL_1_SUBFILTERED",
        #prefix4 = "ALL_1_COLLATED"
        prevFile1 = "Intermediates/COLLATE_{wildcards.location}/" + str(lambda wildcards: get_prev(wildcards.filename)[0]),
        prevFile2 = "Intermediates/COLLATE_{wildcards.location}/" + str(lambda wildcards: get_prev(wildcards.filename)[1])

    run:
        for i in samples:
            shell(f"cp Intermediates/CLEAN/{i}_{wildcards.location}_CLEANED.bed Intermediates/COLLATE_{wildcards.location}/{i}.bed"),
            shell(f"cp Intermediates/CLEAN/{i}_{wildcards.location}_CLEANED.bim Intermediates/COLLATE_{wildcards.location}/{i}.bim"),
            shell(f"cp Intermediates/CLEAN/{i}_{wildcards.location}_CLEANED.fam Intermediates/COLLATE_{wildcards.location}/{i}.fam"),
        datasets = list(config["samples"])
        outputs = datasets
        c = 0
        while len(datasets) > 1:
            datasets = outputs
            outputs = datasets
            hold = ''
            outputName = ''
            c+=1
            if len(datasets) % 2 > 0:
                hold = datasets.pop(-1)
                batch_a = datasets[::2]
                batch_b = datasets[1::2]
                for i,j in zip(batch_a, batch_b):
                    outputName = f"{i}{'-' * c}{j}"
                    shell(f"/apps/plink-1.9/plink --bfile Intermediates/COLLATE_{wildcards.location}/{i} --bmerge Intermediates/COLLATE_{wildcards.location}/{j} --make-bed --keep-allele-order --out Intermediates/COLLATE_{wildcards.location}/{outputName}")
                    del outputs[outputs.index(i)]
                    del outputs[outputs.index(j)]
                    outputs.append(outputName)
                outputs.append(hold)
            else:
                batch_a = datasets[::2]
                batch_b = datasets[1::2]
                for i,j in zip(batch_a, batch_b):
                    outputName = f"{i}{'-' * c}{j}"
                    try:
                        shell(f"/apps/plink-1.9/plink --bfile Intermediates/COLLATE_{wildcards.location}/{i} --bmerge Intermediates/COLLATE_{wildcards.location}/{j} --make-bed --keep-allele-order --out Intermediates/COLLATE_{wildcards.location}/{outputName}")
                    except:
                        print("Tri-Allelic variants found.")
                        if os.path.exists(f"Intermediates/COLLATE_{wildcards.location}/{outputName}-merge.missnp"):
                            print(f"Pulling SNP's from 'Intermediates/COLLATE_{wildcards.location}/{outputName}-merge.missnp'")
                            snps = list()
                            with open(f"Intermediates/COLLATE_{wildcards.location}/{outputName}-merge.missnp", "r") as file:
                                for line in file:
                                    snps.append(line.strip())
                            exclCMD= ",".join(map(str, snps))
                            print(f"Tri-Allelic Variants Identified: {exclCMD}")
                            shell(f"/apps/plink-1.9/plink --bfile Intermediates/COLLATE_{wildcards.location}/{i} --make-bed --keep-allele-order --exclude-snps {exclCMD} --out Intermediates/COLLATE_{wildcards.location}/{i}-FILTERED")
                            shell(f"/apps/plink-1.9/plink --bfile Intermediates/COLLATE_{wildcards.location}/{j} --make-bed --keep-allele-order --exclude-snps {exclCMD} --out Intermediates/COLLATE_{wildcards.location}/{j}-FILTERED")
                            shell(f"/apps/plink-1.9/plink --bfile Intermediates/COLLATE_{wildcards.location}/{i}-FILTERED --bmerge Intermediates/COLLATE_{wildcards.location}/{j}-FILTERED --make-bed --keep-allele-order --out Intermediates/COLLATE_{wildcards.location}/{outputName}")
                    del outputs[outputs.index(str(i))]
                    del outputs[outputs.index(str(j))]
                    outputs.append(outputName)
        fileSet = set([os.path.splitext(file)[0] for file in os.listdir(f"Intermediates/COLLATE_{wildcards.location}/") if all(stringToMatch in file for stringToMatch in config["samples"]) and "merge" not in file])
        for i,j in zip(list(fileSet) * 3, bExtensions):
            shell(f"mv Intermediates/COLLATE_{wildcards.location}/{i}.{j} Intermediates/COLLATE_{wildcards.location}/ALL_{wildcards.location}.{j}")
       
        #module load plink-1.9
        #plink --bfile {params.prevFile1} --bmerge {params.prevFile2} --exclude-snp rs199808813 --make-bed --keep-allele-order --out Intermediates/{wildcards.filename}
       
        #plink --bfile {params.prefix} --bmerge {input.sBed} {input.sBim} {input.sFam} --exclude-snp rs199808813 --make-bed --keep-allele-order --out {params.prefix2}
        #plink --bfile {params.prefix2} --within {input.subPopClusters} --remove-cluster-names ASW ACB  --make-bed --keep-allele-order --out {params.prefix3}
        #plink --bfile {params.prefix3} --make-bed --keep-allele-order --out Intermediates/ALL_{location}

#/*All the *_FILTER processes must output VCF so that I can run through e! Ensembl's Variant Effect Predictor*/
rule ALL_FILTER:
    input:
        #expand("Intermediates/{finalName}.{extension}", extension=["bed", "bim", "fam"], finalName=getFinalName(config["samples"]))
        expand("Intermediates/COLLATE_{{location}}/ALL_{{location}}.{extension}", extension=bExtensions)

    output:
        "Intermediates/FILTER/ALL_{location}_FILTERED.vcf"

    shell:
        """
        module load plink-1.9
        plink --bfile Intermediates/COLLATE_{wildcards.location}/ALL_{wildcards.location} --mind 1 --recode vcf-iid --output-chr chr26 --keep-allele-order --out Intermediates/FILTER/ALL_{wildcards.location}_PRE_SED
        sed -r -e 's/##contig=<ID=chr19,length=[0-9]+>/##contig=<ID=chr19,length=58617616>/' -e 's/##contig=<ID=chr4,length=[0-9]+>/##contig=<ID=chr4,length=190214555>/' Intermediates/FILTER/ALL_{wildcards.location}_PRE_SED.vcf > Intermediates/FILTER/ALL_{wildcards.location}_FILTERED.vcf
        """

rule ALL_ANNOTATE:
    input:
        "Intermediates/FILTER/ALL_{location}_FILTERED.vcf"

    output:
        "Intermediates/ANNOTATE/ALL_{location}_ANNOTATED.vcf"

    shell:
        """
        module load gatk-4.0.12.0
        gatk VariantAnnotator -V {input} -R /apps/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa.gz -D /nlustre/data/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz -O Intermediates/ANNOTATE/ALL_{wildcards.location}_PRE_SED.vcf
        sed -r -e 's/^chr([0-9]{{1,2}})\\t([0-9]+)\\t[0-9]{{1,2}}:[0-9]+[A-Z]{{1}}-[A-Z]{{1}};(rs[0-9]+)/chr\\1\\t\\2\\t\\3/g' Intermediates/ANNOTATE/ALL_{wildcards.location}_PRE_SED.vcf > Intermediates/ANNOTATE/ALL_{wildcards.location}_ANNOTATED.vcf
        """

rule ALL_ANALYZE_SUPER:
    input:
        vcf="Intermediates/ANNOTATE/ALL_{location}_ANNOTATED.vcf",
        popClusters="rawData/superPopCluster"
    
    output:
        expand("Final/SUPER/ALL_{{location}}_SUPER.{extension}", extension=finalExtensions),
        expand("Final/SUPER/{i}/ALL_{{location}}_SUPER_{i}_HV.{extensions}", i=superPop, extensions=["log", "ld"])

    params:
        prefix = 'ALL_{location}_SUPER'

    run:
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --double-id --freq --out Final/SUPER/{params.prefix}"),
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --double-id --within {input.popClusters} --freq --fst --missing --r2 inter-chr dprime --test-mishap --hardy midp --het --ibc --out Final/SUPER/{params.prefix}"),
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --snps-only --double-id --recode HV --out Final/SUPER/ALL_{wildcards.location}_SUPER_HV"),
        shell("module load plink2; plink2 --vcf {input.vcf} --indep-pairwise 50 10 0.1 --double-id --out Final/SUPER/ALL_SUPER_{wildcards.location}"),
        shell("module load plink2; plink2 --vcf {input.vcf} --double-id --mind --extract Final/SUPER/ALL_SUPER_{wildcards.location}.prune.in --pca var-wts scols=sid --out Final/SUPER/ALL_SUPER_{wildcards.location}")
        for i in superPop:
            shell(f"module load plink-1.9; plink --vcf {input.vcf} --double-id --snps-only --keep-allele-order --within {input.popClusters} --keep-cluster-names {i} --r2 inter-chr dprime --recode HV --out Final/SUPER/{i}/{params.prefix}_{i}_HV");
        # Admixture:
        shell("module load plink-1.9; plink --mind --geno --indep-pairwise 50 10 0.1 --vcf {input.vcf} --double-id --keep-allele-order --make-bed --out Final/SUPER/ALL_{wildcards.location}_SUPER"),
        shell("module load admixture-1.3.0; admixture --cv Final/SUPER/ALL_{wildcards.location}_SUPER.bed 5"),
        shell("mv ALL_{wildcards.location}_SUPER.5.* Final/SUPER/")
        
        # /*
        # output:
        # set file("*.bed"), file("*.bim"), file("*.fam") into ALL_1_SUPERPOP
        # file "*.frq" into ALL_1_INITFREQ
        # file "*.frq.strat" into ALL_1_FREQSTRAT

        # plink --vcf {vcfFile} --double-id --within {popClusters} --keep-cluster-names EUR --make-bed --keep-allele-order --out {outPrefix}EUR
        # plink --vcf {vcfFile} --double-id --within {popClusters} --keep-cluster-names EAS --make-bed --keep-allele-order --out {outPrefix}EAS
        # plink --vcf {vcfFile} --double-id --within {popClusters} --keep-cluster-names AMR --make-bed --keep-allele-order --out {outPrefix}AMR
        # plink --vcf {vcfFile} --double-id --within {popClusters} --keep-cluster-names SAS --make-bed --keep-allele-order --out {outPrefix}SAS
        # plink --vcf {vcfFile} --double-id --within {popClusters} --keep-cluster-names AFR --make-bed --keep-allele-order --out {outPrefix}AFR
        # */


rule ALL_ANALYZE_SUB:
    input:
        vcf="Intermediates/FILTER/ALL_{location}_FILTERED.vcf",
        popClusters="rawData/subPopCluster"
        
    output:
        expand("Final/SUB/ALL_{{location}}_SUB.{extension}", extension=finalExtensions),
        expand("Final/SUB/{i}/ALL_{{location}}_SUB_{i}_HV.{extensions}", i=subPop, extensions=["log", "ld"])

    params:
        prefix = 'ALL_{location}_SUB'
  
    run:
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --double-id --freq --out Final/SUB/{params.prefix}"),
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --double-id --within {input.popClusters} --freq --fst --missing --r2 inter-chr dprime --test-mishap --hardy midp --het --ibc --out Final/SUB/{params.prefix}"),
        shell("module load plink-1.9; plink --vcf {input.vcf} --keep-allele-order --snps-only --double-id --recode HV --out Final/SUB/ALL_{wildcards.location}_SUB_HV")
        for i in subPop:
            shell(f"module load plink-1.9; plink --vcf {input.vcf} --double-id --snps-only --keep-allele-order --within {input.popClusters} --keep-cluster-names {i} --r2 inter-chr dprime --recode HV --out Final/SUB/{i}/{params.prefix}_{i}_HV")

    
#// ToDo: Add in Admixture and in-house VEP processes