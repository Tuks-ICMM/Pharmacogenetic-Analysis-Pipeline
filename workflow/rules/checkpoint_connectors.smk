def collect_reportFreqPartitionedPerClusterPerLocation(wildcards):
    checkpoint_output = checkpoints.reportFreqPartitionedPerClusterPerLocation.get(**wildcards).output[0]
    populations = glob_wildcards(join(checkpoint_output, "allele_count.{populations}.acount")).populations
    return expand(outputDir("tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.acount"), cluster=wildcards.cluster, location=wildcards.location, population=populations)

def collect_reportFreqPartitionedPerClusterPerMyLocation(wildcards):
    checkpoint_output = checkpoints.reportFreqPartitionedPerClusterPerLocation.get(**wildcards).output[0]
    return expand(outputDir("tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.acount"), cluster=wildcards.cluster, location=wildcards.location, population=wildcards.population)

def collect_reportHardyWeinburgAllPerLocationOnChrX(wildcards):
    checkpoint_output = checkpoints.reportHardyWeinburgAllPerLocationOnChrX.get(**wildcards).output[0]
    extensions = glob_wildcards(join(checkpoint_output, "%s.{extensions,([\\w\\.]+)}" % wilcards.location)).extensions
    return expand(outputDir("tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.hardy.x"), cluster=wildcards.cluster, location=wildcards.location, population=populations)

def collect_reportHardyWeinburgPartitionedPerClusterPerLocation(wildcards):
    checkpoint_output = checkpoints.reportHardyWeinburgPartitionedPerClusterPerLocation.get(**wildcards).output[0]
    # populations = glob_wildcards(join(checkpoint_output, "hardy_weinberg.{populations}.hardy")).populations
    return expand(outputDir("tmp/reportHardyWeinburgPartitionedPerClusterPerLocation/{cluster}/{location}/hardy_weinberg.{population}.hardy"), cluster=wildcards.cluster, location=wildcards.location, population=wildcards.population)

def collect_reportHardyWeinburgPartitionedPerClusterPerLocationOnChrX(wildcards):
    checkpoint_output = checkpoints.reportHardyWeinburgPartitionedPerClusterPerLocationOnChrX.get(cluster=wildcards.cluster, location=wildcards.location).output[0]
    # populations = glob_wildcards(join(checkpoint_output, "hardy_weinberg.{populations}.hardy.x")).populations
    return expand(outputDir("tmp/reportHardyWeinburgPartitionedPerClusterPerLocationOnChrX/{cluster}/{location}/hardy_weinberg.{population}.hardy.x"), cluster=wildcards.cluster, location=wildcards.location, population=wildcards.population)

def collect_generateSexLinkedTernaryPlotPerClusterPerLocation(wildcards):
    """
    This function is designed to generate a list of ternary plots for sex-linked locations. It is intended for use in the `all` rule.
    """
    files = list()
    for cluster in clusters:
        for location in locations.loc[locations["chromosome"] != 23, "location_name"]:
            for population in samples[cluster].unique():
                files.append(outputDir(f"generateSexLinkedTernaryPlotPerClusterPerLocation/{cluster}/{location}/{population}.jpeg"))
    return files

def collect_generateAutosomalTernaryPlotPerClusterPerLocation(wildcards):
    """
    This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
    """
    files = list()
    for cluster in clusters:
        if cluster in config["fishers-test"]:
            for location in locations.loc[locations["chromosome"] != 23, "location_name"]:
                for population in samples[cluster].unique():
                    files.append(outputDir(f"generateAutosomalTernaryPlotPerClusterPerLocation/{cluster}/{location}/{population}.jpeg"))
    return files


def collect_generateVariantConsequenceBreakdown(wildcards):
    """
    This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
    """
    files = list()
    for cluster in clusters:
        for location in locations["location_name"]:
            for population in samples[cluster].unique():
                files.append(outputDir(f"generateVariantConsequenceBreakdown/{cluster}/{location}/{population}_variant_consequences.jpeg"))
    return files

def collect_generateVariantDistributionByImpact(wildcards):
    """
    This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
    """
    files = list()
    for cluster in clusters:
        for location in locations["location_name"]:
            for population in samples[cluster].unique():
                files.append(outputDir(f"generateVariantDistributionByImpact/{cluster}/{location}/{population}_variant_distribution.jpeg"))
    return files


def collect_analysesFiles(wildcards):
    return [f"tmp/queryVariantEffectPrediction/{wildcards.cluster}/{wildcards.location}/variant_effect_predictions.csv.zst"]