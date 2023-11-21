import os
from io import StringIO
from pathlib import Path

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

        Params.load_from_file(params_file = "tests/pairwise_align.yaml")
        
        # Load parameters
        ref_gene_dir = Params.reference_gene_dir
        ref_gene_file = Params.reference_gene_file
        sequences_very_small = utils.loadSequenceFile(os.path.join(ref_gene_dir, ref_gene_file))
        
        # Build AlignSequence object
        self.sequences = sequences_very_small.copy()
        self.seq_alignment = AlignSequences(self.sequences)
        
        # Get centroid key
        self.centroid = self.seq_alignment.getSequenceCentroid()[0]
        self.centroidSeqs = self.sequences.pop(self.centroid)

        # Get pairwise alignment
        self.aligned = self.seq_alignment.alignSequencesWithPairwise(self.centroid, self.centroidSeqs)
        
    def test_centroidKey(self):
        """
        This test is used to test the centroidKey function.
        """

        print("Begin test_centroidKey...")
        # actual_centroid = self.seq_alignment.getSequenceCentroid()[0]
        
        filename = current_file + "/testFiles/getSequenceCentroid/seq very small"
        with open(filename, "r") as expected_file:
            expected_centroid = expected_file.read()
            assert self.centroid == expected_centroid

    def test_aligned(self):
        """
        This test is used to test the aligned function.
        """

        print("Begin test_aligned...")

        for key in self.aligned.keys():
            filename = current_file + "/testFiles/alignSequence/seq very small/" + key + ".fasta"
            with open(filename, "r") as expected_file:
                expected_file = expected_file.read()
                expected = AlignSequences.fileToDict(current_file + "/testFiles/alignSequence/seq very small/" + key, ".fasta")
                assert self.aligned[key] == expected               
                
    def test_heuristicMSA(self):
        """
        This test is used to test the heuristicMSA function.
        """

        print("Begin test_heuristicMSA...")

        # for alignement in self.alignementSetup:
        starAlignement = self.aligned.starAlignement()
        expected = AlignSequences.fileToDict(current_file + "/testFiles/starAlignement/seq very small", ".fasta")
        assert starAlignement == expected

    def test_windowed(self):
        """
        This test is used to test the windowed function.
        """

        print("Begin test_windowed...")

        for alignement in self.alignementSetup:
            windowed = alignement.windowed

            for key in windowed.keys():
                expected = AlignSequences.fileToDict(current_file + "/testFiles/slidingWindow/seq very small/" + key, ".fasta")
                assert windowed[key] == expected

    def test_msaSet(self):
        """
        This test is used to test the msaSet function.
        """

        print("Begin test_msaSet...")

        for alignement in self.alignementSetup:
            msa = alignement.msa
            
            for key in msa.keys():

                filename = Path(current_file + "/testFiles/makeMSA/seq very small/" + (key + ".fasta"))
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
           
            # Test the createBootstrap function
            genetic_trees = utils.createBoostrap(alignement.msa, p.bootstrapAmount)
            actual_bootstrap = [str(Phylogeny.from_tree(tree)) for tree in list(genetic_trees.values())]

            trees = Phylo.parse("tests/testFiles/createBootstrap/seq very small.xml", "phyloxml")
            expected_bootstrap = [str(tree) for tree in trees]

            for tree in actual_bootstrap:
                assert tree in expected_bootstrap

            # test of the createGeneticList function
            # actual_list, actual_bootstrap_list = utils.createGeneticList(genetic_trees, p.bootstrap_threshold)
            # with open(Path(current_file + "/testFiles/createGeneticList/" + test_case + ".txt"), "r") as f:
            #     expected_list = ast.literal_eval(f.read())
            # assert actual_list == expected_list

            # df = pd.read_csv(p.file_name)
            # climatic_trees = utils.climaticPipeline(df, p.names)
            # utils.filterResults(
            #     climatic_trees, genetic_trees, p.bootstrap_threshold, p.dist_threshold, df, p.reference_gene_filename, p.distance_method
            # )

            # with open(Path(current_file + "/testFiles/writeOutputFiles/" + test_case + ".csv"), "r") as expected_file:
            #     expected_output = [value for value in expected_file.readlines() if value != "\n"]
            # with open("output.csv", "r") as actual_file:
            #     actual_output = [value for value in actual_file.readlines() if value != "\n"]
            # assert len(actual_output) == len(expected_output)
