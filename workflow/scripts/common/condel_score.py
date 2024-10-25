"""
A series of functions needed to calculate a CONDEL score.
"""

from os.path import join
from typing import Tuple

from pandas import read_csv

# %%
#####################################
############ IMPORT DATA ############
#####################################
CONDEL_REFERENCE_DATA = {
    "SIFT": read_csv(
        join("..", "..", "resources", "CONDEL", "sift.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
    "PolyPhen": read_csv(
        join("..", "..", "resources", "CONDEL", "polyphen.data"),
        delimiter="\t",
        names=["Score", "Deleterious", "Normal"],
    ),
}

# %%
################################################
############ SET CONSTANT VARIABLES ############
################################################

# [SET] tests to be used in CONDEL weighted score:
CONDEL_TESTS = ["SIFT", "PolyPhen"]
CONDEL_MINIMUM_CUTOFFS = {"SIFT": 0.15, "PolyPhen": 0.28, "Condel": 0.46}
CONDEL_MAXIMUM_CUTOFFS = {"SIFT": 1, "PolyPhen": 1}


# %%

####################################################
############ DEFINE CONDEL CALCULATIONS ############
####################################################
def condel_weighted_score(sift: int, polyphen: int) -> Tuple[int, str]:
    """Function to calculate a weighted average score accross SIFT and PolyPhenV2 scores.
    Args:
        sift (int): A number indicating the likelyhood a variant affects protein function (Based on sequence homology and protein function)
        polyphen (int): A number indicating impact of a variant on protein structure using straightforward physical comparative methods
    Returns:
        Tuple[int, str]: A weighted score which favours variants with scores further away from each set cutoff (i.e. less ambiguity regarding the prediction) as well as a string, indicating the cutoff verdict
    """

    score = int()

    if sift <= CONDEL_MINIMUM_CUTOFFS["SIFT"]:
        score += (1 - sift / CONDEL_MAXIMUM_CUTOFFS["SIFT"]) * (
            1
            - CONDEL_REFERENCE_DATA["SIFT"]
            .loc[CONDEL_REFERENCE_DATA["SIFT"]["Score"] == sift, "Normal"]
            .values[0]
        )
    else:
        score += (1 - sift / CONDEL_MAXIMUM_CUTOFFS["SIFT"]) * (
            1
            - CONDEL_REFERENCE_DATA["SIFT"]
            .loc[CONDEL_REFERENCE_DATA["SIFT"]["Score"] == sift, "Deleterious"]
            .values[0]
        )

    if polyphen >= CONDEL_MINIMUM_CUTOFFS["PolyPhen"]:
        score += (polyphen / CONDEL_MAXIMUM_CUTOFFS["PolyPhen"]) * (
            1
            - CONDEL_REFERENCE_DATA["PolyPhen"]
            .loc[CONDEL_REFERENCE_DATA["PolyPhen"]["Score"] == polyphen, "Normal"]
            .values[0]
        )
    else:
        score += (polyphen / CONDEL_MAXIMUM_CUTOFFS["PolyPhen"]) * (
            1
            - CONDEL_REFERENCE_DATA["PolyPhen"]
            .loc[CONDEL_REFERENCE_DATA["PolyPhen"]["Score"] == polyphen, "Deleterious"]
            .values[0]
        )

    # This acocunts for number of VEP tests used. Since we are only using 2, we use 2...
    score = score / 2

    if score >= 0.469:
        pred = "deleterious"
    elif score >= 0 and score < 0.469:
        pred = "neutral"

    return (score, pred)
