from os import makedirs
from os.path import exists
from typing import Union
from pandas import DataFrame


def directoryExists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not exists(path):
        makedirs(path)


def search(property: str, rule: str) -> Union[str, int]:
    """Search for a property value defined in the config file, given a property to search for and a rule it should be applied to.

    Args:
        property (str): The name of the property to search for E.g. cores
        rule (str): The name of the rule to search for E.g. VALIDATE

    Returns:
        Union[str, int]: Will return the requested property or error out completely XD.
    """
    return next(
        i[property] for i in config["environment"]["queues"] if rule in i["rules"]
    )

def outputDir(path: str) -> str:
    """This function consults the `config.json` file to determine if a pre-set output directory has been specified. If it has, the provided directory will be used. If not, the current working directory will be used."""
    if "output-dir" in config:
        OUTPUT_DIR_PATH = join(*config["output-dir"])
        return join(OUTPUT_DIR_PATH, path)
    else:
        return join("results", path)

# BEGIN DEFINING RULES:
def vcfValidationWorkflowAdapter(extension: str) -> list:
    """
    An adapter to generate the correct input list from `VCF Validation Pipeline`. This is required as liftover is optional, making the output files variable.
    """
    merge_list = list()
    for reference_genome, groupby_subset in datasets.set_index(["reference_genome", "dataset_name"]).groupby(level=0): # [FOR] all unqiue (dataset_name and reference_genome) column combinations present
        if reference_genome != "GRCh38" and groupby_subset is not None: # [IF] reference genome version
            for dataset_name in groupby_subset.index.get_level_values("dataset_name"): # [FOR] the column in our MultiIndex that contains the dataset_name's in this subset
                # [EACH] add liftover request for the DAG
                merge_list.append(outputDir("tmp/{dataset_name}_liftover{extension}").format(dataset_name=dataset_name, extension=extension))
        else:
            for dataset_name in groupby_subset.index.get_level_values("dataset_name"): # [FOR] the column in our MultiIndex that contains the dataset_name's in this subset
                    # [EACH] add liftover request for the DAG
                    merge_list.append(outputDir("tmp/filter/{dataset_name}_filter{extension}").format(dataset_name=dataset_name, extension=extension))
    return merge_list
