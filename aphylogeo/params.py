import os

import yaml
from yaml.loader import SafeLoader


class Params:
    """
    Class that contains the parameters of the program.
    Loads the parameters from a yaml file or from a dictionary.
    """

    PARAMETER_KEYS = {
        "dist_threshold": 0,
        "window_size": 0,
        "step_size": 0,
        "reference_gene_dir": "",
        "reference_gene_file": "",
        "file_name": "",
        "specimen": "",
        "makeDebugFiles": False,
        "bootstrap_threshold": 0,
        "alignment_method": "0",
        "distance_method": "0",
        "fit_method": "0",
        "tree_type": "0",
        "rate_similarity": 0,
        "method_similarity": "0",
    }

    @classmethod
    def load_from_file(cls, params_file=os.path.join(os.path.dirname(__file__), "params.yaml")):
        """
        Method that loads the parameters from a yaml file.

        args:
            params_file (str): the path to the yaml file

        Raises:
            FileNotFoundError: If the yaml file is not found
        """
        try:
            with open(params_file) as f:
                params = yaml.load(f, Loader=SafeLoader)
                cls.validate_and_set_params(params)
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {params_file}.")

    @classmethod
    def load_default_param(cls):
        cls.validate_and_set_params(cls.PARAMETER_KEYS)

    @classmethod
    def update_from_dict(cls, params_content):
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

        if hasattr(cls, "reference_gene_dir") and hasattr(cls, "reference_gene_file"):
            cls.reference_gene_filepath = os.path.join(cls.reference_gene_dir, cls.reference_gene_file)
        else:
            cls.reference_gene_filepath = None
