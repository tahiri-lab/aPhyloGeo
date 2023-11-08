import os

import yaml
from yaml.loader import SafeLoader


class Params:
    _instance = None

    def __new__(cls, params_file=os.path.join(os.path.dirname(__file__), "params.yaml"), params_content=None):
        if cls._instance is None:
            cls._instance = super(Params, cls).__new__(cls)
            cls._instance.params_file = params_file
            if params_content is None:
                with open(params_file) as f:
                    params = yaml.load(f, Loader=SafeLoader)
            else:
                params = params_content
            cls._instance.__init_params(**params)
        return cls._instance

    def __init_params(self, **params):
        self.bootstrap_threshold = params["bootstrap_threshold"]
        self.dist_threshold = params["dist_threshold"]
        self.window_size = params["window_size"]
        self.step_size = params["step_size"]
        self.data_names = params["data_names"]
        self.reference_gene_file = params["reference_gene_dir"] + "/" + params["reference_gene_file"]
        self.reference_gene_filename = params["reference_gene_file"]

        self.file_name = params["file_name"]
        self.specimen = params["specimen"]
        self.names = params["names"]
        self.makeDebugFiles = params["makeDebugFiles"]
        self.bootstrapAmount = params["bootstrap_amount"]

        self.alignment_method = params["alignment_method"]
        self.distance_method = params["distance_method"]
        self.fit_method = params["fit_method"]
        self.tree_type = params["tree_type"]
        self.rate_similarity = params["rate_similarity"]
        self.method_similarity = params["method_similarity"]
