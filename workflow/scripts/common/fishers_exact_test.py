"""
A series of functions designed to facilitate the calculation of a Fischers Exact Score.
"""

from typing import Tuple

from numpy import nan
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def fishers_exact_test(
    VALUES: dict,
    REFERENCE_POPULATION: str,
    ALTERNATE_POPULATION: str,
    ALTERNATIVE_TEST_NAME: object,
) -> dict:
    """Runs a row-wise Fisher's Exact Test between the two listed populations
    Parameters:
    entry (dict): A dict of DataFrames to work on.
    pop1 (str): A str matching a column name in each dataset corresponding to that populations frequency data.
    compPop (list): A str matching a column name in each dataset corresponding to that populations frequency data.
    """

    contingency = [
        [
            VALUES[f"{REFERENCE_POPULATION}_ac"],
            VALUES[f"{REFERENCE_POPULATION}_tc"] - VALUES[f"{REFERENCE_POPULATION}_ac"],
        ],
        [
            VALUES[f"{ALTERNATE_POPULATION}_ac"],
            VALUES[f"{ALTERNATE_POPULATION}_tc"] - VALUES[f"{ALTERNATE_POPULATION}_ac"],
        ],
    ]

    odds_ratio, p_value = fisher_exact(contingency, alternative=ALTERNATIVE_TEST_NAME)
    return {"OR": odds_ratio, "P": p_value}

    # for fishers_exact_internal_key, fishers_exact_internal_dataset in input_row.items():

    #     fishers_exact_internal_dataset["OR"] = fishers_exact_internal_dataset["Count"][
    #         ["ID", "POS", "REF", "ALT"]
    #     ]
    #     fishers_exact_internal_dataset["P"] = fishers_exact_internal_dataset["Count"][
    #         ["ID", "POS", "REF", "ALT"]
    #     ]
    #     for population in comparison_population:
    #         for key2, direction in {
    #             "L": "less",
    #             "G": "greater",
    #             "T": "two-sided",
    #         }.items():
    #             o_label = f"{reference_population}_{key2}_{population}"
    #             p_label = f"{reference_population}_{key2}_{population}"
    #             fishers_exact_internal_dataset["OR"][o_label] = None
    #             fishers_exact_internal_dataset["P"][p_label] = None
    #             for (
    #                 for_loop_var_index,
    #                 VALUES,
    #             ) in fishers_exact_internal_dataset["Count"].iterrows():
    #                 contingency = [
    #                     [
    #                         VALUES[f"{reference_population}_ac"],
    #                         VALUES[f"{reference_population}_tc"]
    #                         - VALUES[f"{reference_population}_ac"],
    #                     ],
    #                     [
    #                         VALUES[f"{population}_ac"],
    #                         VALUES[f"{population}_tc"]
    #                         - VALUES[f"{population}_ac"],
    #                     ],
    #                 ]
    #                 o_val, p_val = fisher_exact(contingency, alternative=direction)
    #                 fishers_exact_internal_dataset["P"].loc[
    #                     for_loop_var_index, p_label
    #                 ] = p_val
    #                 fishers_exact_internal_dataset["OR"].loc[
    #                     for_loop_var_index, o_label
    #                 ] = o_val
    #             fishers_exact_internal_dataset["P"][p_label] = multipletests(
    #                 fishers_exact_internal_dataset["P"][p_label], method="bonferroni"
    #             )[1]
    # columnsToDrop = list()
    # for pop in populations:
    #     columnsToDrop.append("{}_ac".format(pop))
    #     columnsToDrop.append("{}_tc".format(pop))
    # dataset.drop(columns=columnsToDrop ,inplace=True)
