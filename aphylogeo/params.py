import os

import yaml
from yaml.loader import SafeLoader


class Params:
    """
    Singleton Class that contains the parameters of the program.
    Loads the parameters from a yaml file or from a dictionary.
    """

    _instance = None

    PARAMETER_KEYS = [
        "bootstrap_threshold",
        "dist_threshold",
        "window_size",
        "step_size",
        "data_names",
        "reference_gene_dir",
        "reference_gene_file",
        "file_name",
        "specimen",
        "names",
        "makeDebugFiles",
        "bootstrap_amount",
        "alignment_method",
        "distance_method",
        "fit_method",
        "tree_type",
        "rate_similarity",
        "method_similarity",
    ]

    def __new__(cls):
        """
        Method that creates a singleton of the class.
        """
        if cls._instance is None:
            cls._instance = super(Params, cls).__new__(cls)
        return cls._instance

    @classmethod
    def load_config_from_file(cls, params_file=os.path.join(os.path.dirname(__file__), "params.yaml")):
        """
        Method that loads the parameters from a yaml file.

        args:
            params_file (str): the path to the yaml file
        """
        with open(params_file) as f:
            params = yaml.load(f, Loader=SafeLoader)
            cls.validate_and_set_params(params)

    @classmethod
    def update_config_from_dict(cls, params_content):
        """
        Method that updates the parameters from a dictionary.

        args:
            params_content (dict): the dictionary with the parameters
        """
        cls.validate_and_set_params(params_content)

    @classmethod
    def validate_and_set_params(cls, params_dict):
        """
        Method that validates and sets the parameters.

        args:
            params_dict (dict): the dictionary with the parameters
        """
        for key, value in params_dict.items():
            if key in cls.PARAMETER_KEYS:
                setattr(cls, key, value)
            else:
                raise ValueError(f"Invalid parameter: {key}")
