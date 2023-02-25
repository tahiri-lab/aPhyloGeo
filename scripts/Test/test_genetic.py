from Alignement import AlignSequences
import aPhyloGeo
from Bio import AlignIO
from io import StringIO
from Bio import Phylo
from Bio.Phylo.PhyloXML import Phylogeny
import os
import pandas as pd
from Params import Params
import pytest


current_file = os.path.dirname(__file__)
_params = Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml"))


@pytest.fixture(scope="module")
def climaticTreesSetup():
    return aPhyloGeo.climaticPipeline(_params)

# @pytest.fixture
# def alignementSetup()->list:
#     '''
#     This fixture will be used to create the diffrent alignement objects

#     Returns:
#         list: list of alignement objects
#     '''
#     small = AlignSequences(Params(os.path.join(os.path.dirname(__file__), "params_small.yaml")))
#     very_small = AlignSequences(Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml")))
#     return [very_small, small]

# @pytest.mark.usefixtures('climaticTreesSetup')
# class TestGenetic:

#     def test_genetic_createBootStrap(self):
#         assert True

#     def test_genetic_filterResults(self):
#         assert True




class TestGenetic():

    
    
    def setup_class(self):
        '''
        This fixture will be used to create the diffrent alignement objects

        Returns:
            list: list of alignement objects
        '''
        print("TEST")
        # small = AlignSequences(Params(os.path.join(os.path.dirname(__file__), "params_small.yaml")))
        very_small = AlignSequences(_params)
        self.alignement = [very_small]
        # self.alignementSetup = [very_small, small]


    def test_centroidKey(self):
    
        for alignement in self.alignement:
            
            test_case = alignement.p.reference_gene_filename[0:-6]

            centroid = alignement.centroidKey
            with open(current_file + "\\TestFiles\\GetSequenceCentroid\\" + test_case, 'r') as f:
                centroid_file = f.read()
                assert centroid == centroid_file


    def test_aligned(self):
    
        for alignement in self.alignement:
            
            test_case = alignement.p.reference_gene_filename[0:-6]

            aligned = alignement.aligned
            for key in aligned.keys():
                expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\AlignSequence\\" + test_case + "\\" + key, '.fasta')
                assert aligned[key] == expected


    def test_heuristicMSA(self):
    
        for alignement in self.alignement:
            
            test_case = alignement.p.reference_gene_filename[0:-6]

            starAlignement = alignement.heuristicMSA            
            expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\StarAlignement\\" + test_case, '.fasta')
            assert starAlignement == expected


    def test_windowed(self):

        for alignement in self.alignement:
            
            test_case = alignement.p.reference_gene_filename[0:-6]

            windowed = alignement.windowed
            for key in windowed.keys():
                expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\SlidingWindow\\" + test_case + "\\" + key, '.fasta')
                assert windowed[key] == expected


    def test_msaSet(self):

        for alignement in self.alignement:
            
            test_case = alignement.p.reference_gene_filename[0:-6]

            msa = alignement.msaSet

            for key in msa.keys():
                f=open(current_file + "\\TestFiles\\MakeMSA\\" + test_case + "\\" + key + ".fasta","r")
                data = ""
                noOfLines = 0
                for line in f:
                    data += line
                    noOfLines += 1
                f.close()
                expected = str(AlignIO.read(StringIO(data), "fasta"))
                actual = str(msa[key])

                # for all the lines in expected
                for i in range(noOfLines):
                    if expected[i] not in actual:
                        assert False

    # def test_createBootStrap(self):
    #     for alignement in self.alignement:

    #         test_case = alignement.p.reference_gene_filename[0:-6]

    #         trees = aPhyloGeo.createBoostrap(alignement.msaSet, _params)
    #         # open the file
    #         expected = Phylo.parse(current_file + "\\TestFiles\\CreateBootstrap\\" + test_case + ".xml", "phyloxml")
    #         actual = [str(Phylogeny.from_tree(tree)) for tree in list(trees.values())]
    #         for tree in actual:
    #             print(tree)
            
    #         for tree in expected:
    #             print(str(tree))
    #             if str(tree) not in actual:
    #                 assert False

    def test_filterResults(self):
        assert True


    