#!/usr/bin/env python
"""
A Python script designed to collect the count data from the Plink-2.0 reports. SAves to an Excel file.


[NEEDS] WORKFLOW_RUNTIME

[NEEDED_BY] fishers_exact_bonferonni_corrected.py
[NEEDED_BY] frequency_calculations.py

"""

# pylint: disable=logging-fstring-interpolation
# %%
############ IMPORT DEPENDANCIES ############
#############################################


import logging
from sys import exit

from matplotlib.pyplot import savefig, suptitle
from pandas import merge, read_csv
from upsetplot import plot

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Fatima Barmania",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

# %%


_logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=snakemake.log[0],
    encoding="utf-8",
    filemode="w",
    level=logging.DEBUG,
    format="[%(asctime)s] %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
_logger.info("---LOGGING INTERFACE STARTED---")


# %%
try:

    ############ IMPORT DATA & SET CONSTANTS ############
    #####################################################

    # [SET] standard multi-index Columns
    MULTIINDEX = ["CHROM", "POS", "ID", "REF", "ALT"]

    # [IMPORT] variant count data:
    _logger.info(
        "Importing frequency data for %s | %s from file '%s'",
        snakemake.wildcards.location,
        snakemake.wildcards.cluster,
        snakemake.input.freq_report,
    )
    DATA = read_csv(snakemake.input.freq_report, sep=",")
    _logger.info(
        "Import successful. Available variables are: %s",
        ", ".join(DATA.keys()[:-1]) + f" & {DATA.keys()[-1]}",
    )
    _logger.info(DATA)

    DATA.rename(
        columns={
            population: f"{population}_freq"
            for population in snakemake.params.populations
        },
        inplace=True,
    )
    _logger.debug(
        "Columns renamed successfully. Available keys are: %s",
        ", ".join(DATA.keys()[:-1]) + f" & {DATA.keys()[-1]}",
    )

    # DATA.set_index(MULTIINDEX, inplace=True)
    # _logger.debug("The Multi-index has been created successfully.")

    for population in snakemake.params.populations:
        DATA.loc[(DATA[f"{population}_freq"] < 0.1), population] = True
        DATA.loc[(DATA[f"{population}_freq"] >= 0.1), population] = False
        _logger.debug(
            f"Conversion to boolean type completed for location: {population}"
        )

        DATA.drop(columns=[f"{population}_freq"], inplace=True)
        _logger.debug(f"Old Frequency column for {population} has been dropped.")

    _logger.debug(
        f"Boolean conversion has been completed. Available keys include: {DATA.keys()}"
    )

    _logger.debug(
        f"Attempting boolean multiindexing with available keys: {DATA.keys()} of length {len(DATA.keys())}"
    )
    _logger.debug(DATA)
    DATA.set_index(snakemake.params.populations, inplace=True)
    _logger.debug(
        f"Multi-index has been set successfully for populations: {snakemake.params.populations}."
    )
    _logger.info(DATA)

    plot(
        DATA,
        sort_by="cardinality",
        sort_categories_by="cardinality",
        show_percentages=True,
        subset_size="count",
    )
    _logger.debug("UpPSetPlot has been successfully plotted.")

    suptitle(
        f"{snakemake.wildcards.cluster} | {snakemake.wildcards.location} alleles (<1%) across population intersections",
        fontweight="semibold",
        fontsize="x-large",
    )
    _logger.debug("Figure title has been configured successfully.")

    savefig(
        snakemake.output[0],
    )
    _logger.debug("Plot has been saved successfully.")
# %%

except Exception as E:
    _logger.error(E)
    exit(E)
