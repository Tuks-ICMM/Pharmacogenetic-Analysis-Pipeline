"""
A series of functions designed to facilitate the calculation of a Fischers Exact Score.
"""

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def fishers_exact_test(
    input_row: dict, reference_population: str, comparison_population: list[str]
):
    """Runs a row-wise Fisher's Exact Test between the two listed populations
    Parameters:
    input (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data.
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """
    for fishers_exact_internal_key, fishers_exact_internal_dataset in input_row.items():
        fishers_exact_internal_dataset["OR"] = fishers_exact_internal_dataset["Count"][
            ["ID", "POS", "REF", "ALT"]
        ]
        fishers_exact_internal_dataset["P"] = fishers_exact_internal_dataset["Count"][
            ["ID", "POS", "REF", "ALT"]
        ]
        for population in comparison_population:
            for key2, direction in {
                "L": "less",
                "G": "greater",
                "T": "two-sided",
            }.items():
                o_label = f"{reference_population}_{key2}_{population}"
                p_label = f"{reference_population}_{key2}_{population}"
                fishers_exact_internal_dataset["OR"][o_label] = None
                fishers_exact_internal_dataset["P"][p_label] = None
                for (
                    for_loop_var_index,
                    for_loop_var_row,
                ) in fishers_exact_internal_dataset["Count"].iterrows():
                    contingency = [
                        [
                            for_loop_var_row[f"{reference_population}_ac"],
                            for_loop_var_row[f"{reference_population}_tc"]
                            - for_loop_var_row[f"{reference_population}_ac"],
                        ],
                        [
                            for_loop_var_row[f"{population}_ac"],
                            for_loop_var_row[f"{population}_tc"]
                            - for_loop_var_row[f"{population}_ac"],
                        ],
                    ]
                    o_val, p_val = fisher_exact(contingency, alternative=direction)
                    fishers_exact_internal_dataset["P"].loc[
                        for_loop_var_index, p_label
                    ] = p_val
                    fishers_exact_internal_dataset["OR"].loc[
                        for_loop_var_index, o_label
                    ] = o_val
                fishers_exact_internal_dataset["P"][p_label] = multipletests(
                    fishers_exact_internal_dataset["P"][p_label], method="bonferroni"
                )[1]
        # columnsToDrop = list()
        # for pop in populations:
        #     columnsToDrop.append("{}_ac".format(pop))
        #     columnsToDrop.append("{}_tc".format(pop))
        # dataset.drop(columns=columnsToDrop ,inplace=True)
