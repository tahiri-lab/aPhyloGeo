import ast
import os
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Phylo.PhyloXML import Phylogeny

from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params

current_file = os.path.dirname(__file__)


class TestGenetic:
    def setup_class(self):
        """
        This setup is used to create the diffrent alignement objects, and their
        corresponding params object
        """

        print("Begin setup for test class test_genetic...")

        # params_small = Params(os.path.join(os.path.dirname(__file__), "params_small.yaml"))
        # sequences_small = aPhyloGeo.readFastaFile(params_small.reference_gene_file)
        # small = AlignSequences(params_small.reference_gene_file, params_small.window_size, params_small.step_size,
        #                        params_small.makeDebugFiles, params_small.bootstrapAmount)
        params_very_small = Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml"))
        sequences_very_small = utils.loadSequenceFile(params_very_small.reference_gene_file)
        very_small = AlignSequences(
            sequences_very_small,
            params_very_small.window_size,
            params_very_small.step_size,
            params_very_small.makeDebugFiles,
            params_very_small.bootstrapAmount,
            params_very_small.alignment_method,
            params_very_small.reference_gene_file,
            params_very_small.distance_method,
        )
        very_small.align()
        self.alignementSetup = [very_small]  # , small]
        self.paramSetup = [params_very_small]  # , params_small]

    def test_centroidKey(self):
        """
        This test is used to test the centroidKey function.
        """

        print("Begin test_centroidKey...")

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]
            actual_centroid = alignement.centroidKey
            filename = Path(current_file + "/testFiles/getSequenceCentroid/" + test_case)

            with open(filename, "r") as expected_file:
                expected_centroid = expected_file.read()
                assert actual_centroid == expected_centroid

    def test_aligned(self):
        """
        This test is used to test the aligned function.
        """

        print("Begin test_aligned...")

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]
            aligned = alignement.aligned

            for key in aligned.keys():
                expected = AlignSequences.fileToDict(current_file + "/testFiles/alignSequence/" + test_case + "/" + key, ".fasta")
                assert aligned[key] == expected

    def test_heuristicMSA(self):
        """
        This test is used to test the heuristicMSA function.
        """

        print("Begin test_heuristicMSA...")

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]
            starAlignement = alignement.heuristicMSA
            expected = AlignSequences.fileToDict(current_file + "/testFiles/starAlignement/" + test_case, ".fasta")
            assert starAlignement == expected

    def test_windowed(self):
        """
        This test is used to test the windowed function.
        """

        print("Begin test_windowed...")

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]
            windowed = alignement.windowed

            for key in windowed.keys():
                expected = AlignSequences.fileToDict(current_file + "/testFiles/slidingWindow/" + test_case + "/" + key, ".fasta")
                assert windowed[key] == expected

    def test_msaSet(self):
        """
        This test is used to test the msaSet function.
        """

        print("Begin test_msaSet...")

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]
            msa = alignement.msaSet

            for key in msa.keys():
                filename = Path(current_file + "/testFiles/makeMSA/" + test_case + "/" + (key + ".fasta"))
                f = open(filename, "r")
                data = f.read()
                f.close()

                expected = str(AlignIO.read(StringIO(data), "fasta"))
                actual = str(msa[key])

                for line in expected:
                    assert line in actual

    def test_filterResults(self):
        """
        This test is used to test the filterResults function.
        """

        for alignement, p in zip(self.alignementSetup, self.paramSetup):
            test_case = p.reference_gene_filename[0:-6]

            # Test the createBootstrap function
            genetic_trees = utils.createBoostrap(alignement.msaSet, p.bootstrapAmount)
            actual_bootstrap = [str(Phylogeny.from_tree(tree)) for tree in list(genetic_trees.values())]
            actual_bootstrap = [(tree.splitlines()).sort() for tree in actual_bootstrap]

            expected_bootstrap = [str(tree) for tree in Phylo.parse(current_file + "/testFiles/createBootstrap/" + test_case + ".xml", "phyloxml")]
            expected_bootstrap = [(tree.splitlines()).sort() for tree in expected_bootstrap]

            for tree in actual_bootstrap:
                assert tree in expected_bootstrap

            # test of the createGeneticList function
            actual_list, actual_bootstrap_list = utils.createGeneticList(genetic_trees, p.bootstrap_threshold)
            with open(Path(current_file + "/testFiles/createGeneticList/" + test_case + ".txt"), "r") as f:
                expected_list = ast.literal_eval(f.read())
            assert actual_list == expected_list

            df = pd.read_csv(p.file_name)
            climatic_trees = utils.climaticPipeline(df, p.names)
            utils.filterResults(
                climatic_trees, genetic_trees, p.bootstrap_threshold, p.dist_threshold, df, p.reference_gene_filename, p.distance_method
            )

            with open(Path(current_file + "/testFiles/writeOutputFiles/" + test_case + ".csv"), "r") as expected_file:
                expected_output = [value for value in expected_file.readlines() if value != "\n"]
            with open("output.csv", "r") as actual_file:
                actual_output = [value for value in actual_file.readlines() if value != "\n"]
            assert len(actual_output) == len(expected_output)
