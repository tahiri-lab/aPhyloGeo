import json

from Bio import Phylo


class GeneticTrees:
    """
    Class that contains the data of Window.
    """

    def __init__(self, trees_dict: dict, format: str = "newick"):
        """
        init method
        args:
            trees_dict (dict): {str, biopython tree}
            format (str): the format of the trees - default: newick
        """
        self.trees = trees_dict
        self.format = format

    def get_trees_str(self):
        """
        Get the trees in string format.

        Return:
            trees (dict): {str, str}
        """
        trees = {}
        for key, value in self.trees.items():
            trees[key] = value.format(self.format)

        return trees

    def save_trees_to_json(self, file_name: str, format: str = "newick"):
        """
        Save the trees in json format.

        args:
            file_name (str): the name of the file
            format (str): the format of the trees - default: newick
        """
        genetic_trees = {}
        for key, value in self.trees.items():
            genetic_trees[key] = value.format(format)
        with open(file_name, "w") as f:
            json.dump(genetic_trees, f)

    @classmethod
    def load_trees_from_json(cls, file_name: str, format: str = "newick"):
        """
        Load the trees from json format.

        args:
            file_name (str): the name of the file
            format (str): the format of the trees - default: newick
        """
        with open(file_name, "r") as f:
            genetic_trees = json.load(f)
        trees = {}
        for key, value in genetic_trees.items():
            trees[key] = Phylo.read(value, format)
        return cls(trees)
