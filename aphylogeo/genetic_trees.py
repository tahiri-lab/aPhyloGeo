import json
import io

from Bio import Phylo


class GeneticTrees:
    """
    Class that contains the genetic trees. It can save and load the trees from file or json format.
    The trees are stored in a dictionary {str(window), biopython tree}
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

    def get_trees_str(self) -> dict:
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
    def load_trees_from_file(cls, file_name: str, format: str = "newick") -> dict:
        """
        Load the trees from json format.

        args:
            file_name (str): the name of the file
            format (str): the format of the trees - default: newick
                Can be called with None as format and Biopython will identitfy the format.

        return:
            GeneticTrees: dict of windows, with the trees in biopython format.
        """
        with open(file_name, "r") as f:
            str_trees = f.read()
        trees = cls.load_trees_from_json(str_trees)
        return trees

    @classmethod
    def load_trees_from_json(cls, trees: str, format: str = "newick") -> dict:
        """
        Load the trees from json format.

        args:
            trees (str): the trees in json format
            format (str): the format of the trees - default: newick
                Can be called with None as format and Biopython will identitfy the format.

        return:
            GeneticTrees: dict of windows, with the trees in biopython format.
        """

        genetic_trees = json.loads(trees)
        trees = {}
        for key, value in genetic_trees.items():
            with io.StringIO(value.strip()) as file:
                trees[key] = Phylo.read(file, format=format)
        return cls(trees)

    @classmethod
    def testtrees(cls, file_name: str, format: str = "newick") -> dict:
        """
        Load the trees from json format.

        args:
            file_name (str): the name of the file
            format (str): the format of the trees - default: newick
                Can be called with None as format and Biopython will identitfy the format.
        """
        with open(file_name, "r") as f:
            str_trees = f.read()
        trees = cls.load_trees_from_json(str_trees)
        return cls(trees)
