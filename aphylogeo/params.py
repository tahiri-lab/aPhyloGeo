import yaml
from yaml.loader import SafeLoader
import os


class Params:

    def __init__(self, params_file=os.path.join(os.path.dirname(__file__), "params.yaml")):

        # We open the params.yaml file and put it in the params variable
        # We use the absolute path to the file so that we can run the script from anywhere
        with open(params_file) as f:
            params = yaml.load(f, Loader=SafeLoader)

        self.bootstrap_threshold = params["bootstrap_threshold"]
        self.dist_threshold = params["dist_threshold"]
        self.window_size = params["window_size"]
        self.step_size = params["step_size"]
        self.data_names = params["data_names"]
        self.reference_gene_file = params["reference_gene_dir"] + '/' + params["reference_gene_file"]
        self.reference_gene_filename = params["reference_gene_file"]

        self.file_name = params["file_name"]
        self.specimen = params["specimen"]
        self.names = params["names"]
        self.makeDebugFiles = params["makeDebugFiles"]
        self.bootstrapAmount = params["bootstrap_amount"]

        self.alignment_method = params["alignment_method"]
        self.distance_method = params["distance_method"]
        # self.fasta_path = params["fasta_path"]
