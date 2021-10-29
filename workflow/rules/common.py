from os.path import exists
from os import makedirs
from typing import Union


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
