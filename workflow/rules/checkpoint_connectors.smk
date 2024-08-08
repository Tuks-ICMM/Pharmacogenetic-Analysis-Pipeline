def collect_reportFreqPartitionedPerClusterPerLocation(wildcards):
    checkpoint_output = checkpoints.reportFreqPartitionedPerClusterPerLocation.get(
        **wildcards
    ).output[0]
    populations = glob_wildcards(
        join(checkpoint_output, "allele_count.{populations}.acount")
    ).populations
    return expand(
        outputDir(
            "tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.acount"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


def collect_reportLinkageDisequilibriumPartitionedPerClusterPerLocation(wildcards):
    FILES_TO_RETURN = list()
    for location in locations["location_name"].tolist():
        for cluster in clusters:
            checkpoint_output = checkpoints.reportLinkageDisequilibriumPartitionedPerClusterPerLocation.get(
                location=location, cluster=cluster
            ).output["linkage_report"]

            FILES_TO_RETURN.append(checkpoint_output)
    return FILES_TO_RETURN


def collect_reportFreqPartitionedPerClusterPerMyLocation(wildcards):
    checkpoint_output = checkpoints.reportFreqPartitionedPerClusterPerLocation.get(
        **wildcards
    ).output[0]
    return expand(
        outputDir(
            "tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.acount"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=wildcards.population,
    )


def collect_reportHardyWeinbergAllPerLocationOnChrX(wildcards):
    checkpoint_output = checkpoints.reportHardyWeinburgAllPerLocationOnChrX.get(
        **wildcards
    ).output[0]
    extensions = glob_wildcards(
        join(checkpoint_output, "%s.{extensions,([\\w\\.]+)}" % wilcards.location)
    ).extensions
    return expand(
        outputDir(
            "tmp/reportFreqPartitionedPerClusterPerLocation/{cluster}/{location}/allele_count.{population}.hardy.x"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


def collect_reportAutosomalHardyWeinbergPartitionedPerClusterPerLocation(wildcards):
    checkpoint_output = (
        checkpoints.reportAutosomalHardyWeinbergPartitionedPerClusterPerLocation.get(
            **wildcards
        ).output[0]
    )
    populations = glob_wildcards(
        join(checkpoint_output, "hardy_weinberg.{populations}.hardy")
    ).populations
    return expand(
        outputDir(
            "tmp/reportAutosomalHardyWeinbergPartitionedPerClusterPerLocation/{cluster}/{location}/hardy_weinberg.{population}.hardy"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


def collect_reportAutosomalHardyWeinbergPartitionedPerClusterPerLocationOnChrX(
    wildcards,
):
    checkpoint_output = checkpoints.reportAutosomalHardyWeinbergPartitionedPerClusterPerLocationOnChrX.get(
        cluster=wildcards.cluster, location=wildcards.location
    ).output[
        0
    ]
    # populations = glob_wildcards(join(checkpoint_output, "hardy_weinberg.{populations}.hardy.x")).populations
    return expand(
        outputDir(
            "tmp/reportAutosomalHardyWeinbergPartitionedPerClusterPerLocationOnChrX/{cluster}/{location}/hardy_weinberg.{population}.hardy.x"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=wildcards.population,
    )


def collect_reportFixationIndexPerClusterPerLocation(wildcards):
    FILES_TO_RETURN = list()

    for location in locations.loc[
        (locations["chromosome"] <= 22) & (locations["chromosome"] >= 1),
        "location_name",
    ]:
        for cluster in clusters:
            FILES_TO_RETURN.append(
                checkpoints.reportFixationIndexPerClusterPerLocation.get(
                    cluster=cluster, location=location
                ).output["fixation_report"]
            )
    return FILES_TO_RETURN


def collect_reportSampleMissingnessPerClusterPerLocation(wildcards):
    FILES_TO_RETURN = list()
    for location in locations["location_name"]:
        for cluster in clusters:
            checkpoint_output = checkpoints.reportMissingnessPerClusterPerLocation.get(
                cluster=cluster, location=location
            ).output[0]
            populations = glob_wildcards(
                join(
                    checkpoint_output,
                    "{cluster}_{location}_missingness.{populations}.smiss.zst",
                )
            ).populations
            FILES_TO_RETURN.extend(
                expand(
                    outputDir(
                        "tmp/reportMissingnessPerClusterPerLocation/{cluster}/{location}/{cluster}_{location}_missingness.{population}.smiss.zst"
                    ),
                    cluster=cluster,
                    location=location,
                    population=populations,
                )
            )
    return FILES_TO_RETURN



def collect_reportVariantMissingnessPerClusterPerLocation(wildcards):
    checkpoint_output = checkpoints.reportMissingnessPerClusterPerLocation.get(
        **wildcards
    ).output[0]
    populations = glob_wildcards(
        join(
            checkpoint_output,
            "{cluster}_{location}_missingness.{populations}.vmiss.zst",
        )
    ).populations
    return expand(
        outputDir(
            "tmp/reportMissingnessPerClusterPerLocation/{cluster}/{location}/{cluster}_{location}_missingness.{population}.vmiss.zst"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


# def collect_generateSexLinkedTernaryPlotPerClusterPerLocation(wildcards):
#     """
#     This function is designed to generate a list of ternary plots for sex-linked locations. It is intended for use in the `all` rule.
#     """
#     files = list()
#     for cluster in clusters:
#         for location in locations.loc[locations["chromosome"] != 23, "location_name"]:
#             for population in samples[cluster].unique():
#                 files.append(
#                     outputDir(
#                         f"generateSexLinkedTernaryPlotPerClusterPerLocation/{cluster}/{location}/{population}.jpeg"
#                     )
#                 )
#     return files


# def collect_generateAutosomalTernaryPlotPerClusterPerLocation(wildcards):
#     """
#     This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
#     """
#     files = list()
#     for cluster in clusters:
#         if cluster in config["fishers-test"]:
#             for location in locations.loc[
#                 locations["chromosome"] != 23, "location_name"
#             ]:
#                 for population in samples[cluster].unique():
#                     files.append(
#                         outputDir(
#                             f"generateAutosomalTernaryPlotPerClusterPerLocation/{cluster}/{location}/{population}.jpeg"
#                         )
#                     )
#     return files


# def collect_generateVariantConsequenceBreakdown(wildcards):
#     """
#     This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
#     """
#     files = list()
#     for cluster in clusters:
#         for location in locations["location_name"]:
#             for population in samples[cluster].unique():
#                 files.append(
#                     outputDir(
#                         f"generateVariantConsequenceBreakdown/{cluster}/{location}/{population}_variant_consequences.jpeg"
#                     )
#                 )
#     return files


# def collect_generateVariantDistributionByImpact(wildcards):
#     """
#     This function is designed to generate a list of ternary plot files. It is intended for use in the `all` rule.
#     """
#     files = list()
#     for cluster in clusters:
#         for location in locations["location_name"]:
#             for population in samples[cluster].unique():
#                 files.append(
#                     outputDir(
#                         f"generateVariantDistributionByImpact/{cluster}/{location}/{population}_variant_distribution.jpeg"
#                     )
#                 )
#     return files


# def collect_calculateLinkageDisequilibriumPartitionedPerClusterPerLocation(wildcards):
#     """
#     This function is designed to colelct all cluster-level outputs from the `calculateLinkageDisequilibriumPartitionedPerClusterPerLocation` rule.
#     """
#     files = list()
#     for cluster in clusters:
#         for location in locations["location_name"]:
#             for population in samples[cluster].unique():
#                 files.append(
#                     outputDir(
#                         f"generateVariantDistributionByImpact/{cluster}/{location}/{population}_variant_distribution.jpeg"
#                     )
#                 )
#     return files


def collect_filesToConsolidate(wildcards):
    base = [
        outputDir(
            f"tmp/compileVariantEffectPrediction/{wildcards.cluster}/{wildcards.location}/cleaned_variant_effect_predictions.csv.zst"
        ),
        outputDir(
            f"tmp/collectFreqPartitionedPerClusterPerLocation/{wildcards.cluster}_{wildcards.location}.csv.zst"
        ),
        outputDir(
            f"tmp/collectCountPartitionedPerClusterPerLocation/{wildcards.cluster}_{wildcards.location}.csv.zst"
        ),
        outputDir(
            f"tmp/collectAutosomalHardyWeinbergPartitionedPerClusterPerLocation/{wildcards.cluster}/{wildcards.location}/{wildcards.cluster}_{wildcards.location}_hardy_weinberg.csv.zst"
        ),
        outputDir(f"tmp/collectVariantMissingnessPerClusterPerLocation/{wildcards.cluster}/{wildcards.location}/{wildcards.cluster}_{wildcards.location}_vmiss.csv.zst")
    ]

    if wildcards.cluster in config["fishers-test"]:
        base.append(
            outputDir(
                f"tmp/calculateFichersExactTestWithCorrection/{wildcards.cluster}/{wildcards.location}/fishers_exact.csv.zst"
            )
        )
    return base
