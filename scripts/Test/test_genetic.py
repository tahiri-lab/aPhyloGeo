from Alignement import AlignSequences
import aPhyloGeo
import ast
from Bio import AlignIO
from io import StringIO
from Bio import Phylo
from Bio.Phylo.PhyloXML import Phylogeny
import os
from Params import Params
from pathlib import Path


current_file = os.path.dirname(__file__)


class TestGenetic():

    def setup_class(self):
        '''
        This setup is used to create the diffrent alignement objects, and their
        corresponding params object
        '''

        print("Begin setup for test class test_genetic...")

        params_small = Params(os.path.join(os.path.dirname(__file__), "params_small.yaml"))
        params_very_small = Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml"))
        small = AlignSequences(params_small)
        very_small = AlignSequences(params_very_small)

        self.alignementSetup = [very_small, small]
        self.paramSetup = [params_very_small, params_small]

    def test_centroidKey(self):
        ''' 
        This test is used to test the centroidKey function.
        '''

        print("Begin test_centroidKey...")
    
        for alignement in self.alignementSetup:
            
            test_case = alignement.p.reference_gene_filename[0:-6]
            actual_centroid = alignement.centroidKey
            filename = Path(current_file + "/TestFiles/GetSequenceCentroid/" + test_case)

            with open(filename, 'r') as expected_file:
                expected_centroid = expected_file.read()
                assert actual_centroid == expected_centroid

    def test_aligned(self):
        '''
        This test is used to test the aligned function.
        '''

        print("Begin test_aligned...")
    
        for alignement in self.alignementSetup:
            
            test_case = alignement.p.reference_gene_filename[0:-6]
            aligned = alignement.aligned

            for key in aligned.keys():
                expected = AlignSequences.fileToDict(current_file + "/TestFiles/AlignSequence/" + test_case + "/" + key, '.fasta')
                assert aligned[key] == expected

    def test_heuristicMSA(self):
        '''
        This test is used to test the heuristicMSA function.
        '''

        print("Begin test_heuristicMSA...")
    
        for alignement in self.alignementSetup:
            
            test_case = alignement.p.reference_gene_filename[0:-6]
            starAlignement = alignement.heuristicMSA            
            expected = AlignSequences.fileToDict(current_file + "/TestFiles/StarAlignement/" + test_case, '.fasta')
            assert starAlignement == expected

    def test_windowed(self):
        '''
        This test is used to test the windowed function.
        '''

        print("Begin test_windowed...")

        for alignement in self.alignementSetup:
            
            test_case = alignement.p.reference_gene_filename[0:-6]
            windowed = alignement.windowed

            for key in windowed.keys():
                expected = AlignSequences.fileToDict(current_file + "/TestFiles/SlidingWindow/" + test_case + "/" + key, '.fasta')
                assert windowed[key] == expected

    def test_msaSet(self):
        '''
        This test is used to test the msaSet function.
        '''

        print("Begin test_msaSet...")

        for alignement in self.alignementSetup:
            
            test_case = alignement.p.reference_gene_filename[0:-6]
            msa = alignement.msaSet
            
            for key in msa.keys():

                filename = Path(current_file + "/TestFiles/MakeMSA/" + test_case + "/" + (key + ".fasta"))
                f = open(filename, "r")
                data = f.read()
                f.close()

                expected = str(AlignIO.read(StringIO(data), "fasta"))
                actual = str(msa[key])

                for line in expected:
                    assert line in actual

    def test_filterResults(self):
        '''
        This test is used to test the filterResults function.
        '''

        for alignement, params in zip(self.alignementSetup, self.paramSetup):

            test_case = alignement.p.reference_gene_filename[0:-6]

            # Test the createBootstrap function
            genetic_trees = aPhyloGeo.createBoostrap(alignement.msaSet, params)
            actual_bootstrap = [str(Phylogeny.from_tree(tree)) for tree in list(genetic_trees.values())]
            actual_bootstrap = [(tree.splitlines()).sort() for tree in actual_bootstrap]
            
            expected_bootstrap = [str(tree) for tree in Phylo.parse(current_file + "/TestFiles/CreateBootstrap/" + test_case + ".xml", "phyloxml")]
            expected_bootstrap = [(tree.splitlines()).sort() for tree in expected_bootstrap]
            
            for tree in actual_bootstrap:
                assert tree in expected_bootstrap

            # test of the createGeneticList function
            actual_list, actual_bootstrap_list = aPhyloGeo.createGeneticList(genetic_trees, params)
            with open(Path(current_file + "/TestFiles/CreateGeneticList/" + test_case + ".txt"), 'r') as f:
                expected_list = ast.literal_eval(f.read())
            assert actual_list == expected_list

            climatic_trees = aPhyloGeo.climaticPipeline(params)
            aPhyloGeo.filterResults(climatic_trees, genetic_trees, params)

            # TO BE CONFIRMED
            # test of the writeOutputFiles function
            with open(Path(current_file + "/TestFiles/WriteOutputFiles/" + test_case + ".csv"), 'r') as expected_file:
                expected_output = expected_file.readlines().remove("\n")
            with open("output.csv", 'r') as actual_file:
                actual_output = actual_file.readlines().remove("\n")
            assert len(actual_output) == len(expected_output)
