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
