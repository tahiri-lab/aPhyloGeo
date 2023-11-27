import os
from io import StringIO
from pathlib import Path
import glob

import Bio.SeqIO
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

        # Get Heuristic MSA
        self.heuristicMSA = self.seq_alignment.starAlignement(self.centroid, self.aligned)

        # Get windowed alignment
        self.windowed = self.seq_alignment.slidingWindow(self.heuristicMSA)

        # MakeMSA
        self.msa = self.seq_alignment.makeMSA(self.windowed)

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
        expected = AlignSequences.fileToDict(current_file + "/testFiles/starAlignement/seq very small", ".fasta")
        assert self.heuristicMSA == expected

    def test_windowed(self):
        """
        This test is used to test the windowed function.
        """

        print("Begin test_windowed...")

        for key in self.windowed.keys():
            expected = AlignSequences.fileToDict(current_file + "/testFiles/slidingWindow/seq very small/" + key, ".fasta")
            assert self.windowed[key] == expected

    def test_msaSet(self):
        """
        This test is used to test the msaSet function.
        """

        print("Begin test_msaSet...")
            
        for key in self.msa.keys():

            filename = Path(current_file + "/testFiles/makeMSA/seq very small/" + (key + ".fasta"))
            f = open(filename, "r")
            data = f.read()
            f.close()

            expected = str(AlignIO.read(StringIO(data), "fasta"))
            actual = str(self.msa[key])

            for line in expected:
                assert line in actual

    def test_filterResults(self):
        """
        This test is used to test the filterResults function.
        """

        # Test the createBootstrap function
        genetic_trees = utils.createBoostrap(self.msa, Params.bootstrap_amount)
        actual_bootstrap = [str(Phylogeny.from_tree(tree)) for tree in list(genetic_trees.values())]

        trees = Phylo.parse("tests/testFiles/createBootstrap/seq very small.xml", "phyloxml")
        expected_bootstrap = [str(tree) for tree in trees]

        for tree in actual_bootstrap:
            assert tree in expected_bootstrap

    def test_clustal(self):
        """
        This test is used to test the clustalAlign function.
        """
        Params.load_from_file(params_file = "tests/clustal_align.yaml")

        # load parameters
        ref_gene_dir = Params.reference_gene_dir
        ref_gene_file = Params.reference_gene_file
        sequences_very_small = utils.loadSequenceFile(os.path.join(ref_gene_dir, ref_gene_file))
        
        # Build AlignSequence object
        self.sequences = sequences_very_small.copy()
        self.seq_alignment = AlignSequences(self.sequences)

        # Call clustal Alignment
        clustal_alignment = self.seq_alignment.clustalAlign()
        [os.remove(file) for file in glob.glob("bin/tmp/*.fasta")]

        # Load expected alignment
        fasta_out = "tests/testFiles/clustalAlign/clustal_alignment.fasta"
        records = Bio.SeqIO.parse(fasta_out, "fasta")
        expected = {rec.id: str(rec.seq) for rec in records}
        
        # Assert
        assert clustal_alignment == expected

    def test_muscle(self):
        """
        This test is used to test the muscleAlign function.
        """
        Params.load_from_file(params_file = "tests/muscle_align.yaml")

        # load parameters
        ref_gene_dir = Params.reference_gene_dir
        ref_gene_file = Params.reference_gene_file
        sequences_very_small = utils.loadSequenceFile(os.path.join(ref_gene_dir, ref_gene_file))
        
        # Build AlignSequence object
        self.sequences = sequences_very_small.copy()
        self.seq_alignment = AlignSequences(self.sequences)

        # Call muscle Alignment
        muscle_alignment = self.seq_alignment.muscleAlign()
        [os.remove(file) for file in glob.glob("bin/tmp/*.fasta")]

        # Load expected alignment
        fasta_out = "tests/testFiles/muscleAlign/muscle_alignment.fasta"
        records = Bio.SeqIO.parse(fasta_out, "fasta")
        expected = {rec.id: str(rec.seq) for rec in records}

        # Assert
        assert muscle_alignment == expected

    def test_mafft(self):
        """
        This test is used to test the mafftAlign function.
        """
        Params.load_from_file(params_file = "tests/mafft_align.yaml")

        # load parameters
        ref_gene_dir = Params.reference_gene_dir
        ref_gene_file = Params.reference_gene_file
        sequences_very_small = utils.loadSequenceFile(os.path.join(ref_gene_dir, ref_gene_file))
        
        # Build AlignSequence object
        self.sequences = sequences_very_small.copy()
        self.seq_alignment = AlignSequences(self.sequences)

        # Call muscle Alignment
        muscle_alignment = self.seq_alignment.mafftAlign()
        
        # Load expected alignment
        fasta_out = "tests/testFiles/mafftAlign/mafft_alignment.fasta"
        records = Bio.SeqIO.parse(fasta_out, "fasta")
        expected = {rec.id: str(rec.seq) for rec in records}

        # Assert
        assert muscle_alignment == expected