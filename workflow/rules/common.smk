__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "rob@spot.colorado.edu"
__status__ = "Beta"


def directoryExists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not os.path.exists(path):
        os.makedirs(path)
