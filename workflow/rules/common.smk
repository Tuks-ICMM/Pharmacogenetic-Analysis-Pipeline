def directoryExists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not os.path.exists(path):
                    os.makedirs(path)