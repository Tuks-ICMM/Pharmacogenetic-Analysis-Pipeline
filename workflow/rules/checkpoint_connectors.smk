
def collect_report_freq_partitioned_per_cluster(wildcards):
    checkpoint_output = checkpoints.report_freq_partitioned_per_cluster.get(
        **wildcards
    ).output["files"]
    glob_match = glob_wildcards(
        join(checkpoint_output, "allele_count.{populations}.acount")
    )
    return expand(
        outputDir(
            "tmp/{cluster}/{location}/freqByCluster/allele_count.{population}.acount"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=glob_match.populations,
    )


def collect_calculate_linkage_disequilibrium_per_cluster(wildcards) -> list[str]:
    LD_output = list()
    for location in locations["location_name"].unique().tolist():
        for cluster in clusters:
            LD_output.append(
                directory(outputDir(f"linkage_disequilibrium/{cluster}/{location}/"))
            )
    return LD_output


# def collect_reportFreqPartitionedPerClusterPerMyLocation(wildcards):
#     checkpoint_output = checkpoints.report_freq_partitioned_per_cluster.get(
#         **wildcards
#     ).output[0]
#     return expand(
#         outputDir(
#             "tmp/report_freq_partitioned_per_cluster/{cluster}/{location}/allele_count.{population}.acount"
#         ),
#         cluster=wildcards.cluster,
#         location=wildcards.location,
#         population=wildcards.population,
#     )


def collect_sex_linked_hardy_weinberg_per_cluster(wildcards):
    checkpoint_output = checkpoints.reportHardyWeinbergAllPerLocationOnChrX.get(
        **wildcards
    ).output[0]
    extensions = glob_wildcards(
        join(checkpoint_output, "%s.{extensions,([\\w\\.]+)}" % wilcards.location)
    ).extensions
    return expand(
        outputDir(
            "tmp/report_freq_partitioned_per_cluster/{cluster}/{location}/allele_count.{population}.hardy.x"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


def collect_autosomal_hardy_weinberg_per_cluster(wildcards):
    checkpoint_output = checkpoints.report_hardy_weinberg_per_cluster.get(**wildcards).output["directory"]
    populations = glob_wildcards(
        join(checkpoint_output, "hardy_weinberg.{populations}.hardy")
    ).populations
    return expand(
        outputDir(
            "tmp/{cluster}/{location}/hardy_weinberg_per_cluster/hardy_weinberg.{population}.hardy"
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


def collect_report_fixation_index_per_cluster(wildcards):
    FILES_TO_RETURN = list()

    for location in locations.loc[
        (locations["chromosome"] <= 22) & (locations["chromosome"] >= 1),
        "location_name",
    ]:
        for cluster in clusters:
            # FILES_TO_RETURN.append(
            #     checkpoints.report_fixation_index_per_cluster.get(
            #         cluster=cluster, location=location
            #     ).output["fixation_report"]
            # )
            FILES_TO_RETURN.append(
                directory(outputDir(f"fixation_index/{cluster}/{location}/"))
            )
    return FILES_TO_RETURN


def collect_reportSampleMissingnessPerClusterPerLocation(wildcards):
    FILES_TO_RETURN = list()
    for location in locations["location_name"]:
        for cluster in clusters:
            checkpoint_output = checkpoints.report_missingness_per_cluster.get(
                cluster=cluster, location=location
            ).output[0]
            populations = glob_wildcards(
                join(
                    checkpoint_output,
                    "{cluster}_{location}_missingness.{populations}.smiss",
                )
            ).populations
            FILES_TO_RETURN.extend(
                expand(
                    outputDir(
                        "tmp/report_missingness_per_cluster/{cluster}/{location}/{cluster}_{location}_missingness.{population}.smiss"
                    ),
                    cluster=cluster,
                    location=location,
                    population=populations,
                )
            )
    return FILES_TO_RETURN


def collect_report_missingness_per_cluster(wildcards):
    checkpoint_output = checkpoints.report_missingness_per_cluster.get(
        **wildcards
    ).output[0]
    populations = glob_wildcards(
        join(
            checkpoint_output,
            "{cluster}_{location}_missingness.{populations}.vmiss",
        )
    ).populations
    return expand(
        outputDir(
            "tmp/report_missingness_per_cluster/{cluster}/{location}/{cluster}_{location}_missingness.{population}.vmiss"
        ),
        cluster=wildcards.cluster,
        location=wildcards.location,
        population=populations,
    )


def collect_reports_to_consolidate(wildcards):
    base = [
        outputDir(
            f"tmp/{wildcards.cluster}/{wildcards.location}/cleaned_variant_effect_predictions.csv"
        ),
        outputDir(
            f"tmp/{wildcards.cluster}/{wildcards.location}/variant_frequency.csv"
        ),
        outputDir(
            f"tmp/{wildcards.cluster}/{wildcards.location}/variant_count.csv"
        ),
        outputDir(
            f"tmp/{wildcards.cluster}/{wildcards.location}/autosomal_hardy_weinberg.csv"
        ),
        outputDir(f"tmp/{wildcards.cluster}/{wildcards.location}/missingness.csv"),
    ]

    if wildcards.cluster in config["fishers-test"]:
        base.append(
            outputDir(
                f"tmp/{wildcards.cluster}/{wildcards.location}/fishers_exact_with_corrections.csv"
            )
        )
    return base
