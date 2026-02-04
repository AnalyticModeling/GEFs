import os
def get_path():
    # Returns the absolute path to the matlab_wrappers folder
    return os.path.join(os.path.dirname(__file__), 'matlab_wrappers')