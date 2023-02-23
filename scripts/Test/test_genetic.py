from Alignement import AlignSequences
import aPhyloGeo
from Bio import AlignIO
from io import StringIO
import os
import pandas as pd
from Params import Params
import pytest


current_file = os.path.dirname(__file__)

# @pytest.fixture(scope="module")
# def climaticTreesSetup():
#     return aPhyloGeo.climaticPipeline()

@pytest.fixture
def alignementSetup()->list:
    '''
    This fixture will be used to create the diffrent alignement objects

    Returns:
        list: list of alignement objects
    '''
    small = AlignSequences(Params(os.path.join(os.path.dirname(__file__), "params_small.yaml")))
    very_small = AlignSequences(Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml")))
    return [very_small, small]

# @pytest.mark.usefixtures('climaticTreesSetup')
# class TestGenetic:

#     def test_genetic_alignement(self, climaticTreesSetup):
#         assert True

#     def test_genetic_createBootStrap(self):
#         assert True

#     def test_genetic_filterResults(self):
#         assert True


@pytest.mark.usefixtures('alignementSetup')
def test_alignement(alignementSetup):
    
    for alignement in alignementSetup:
        
        test_case = alignement.p.reference_gene_filename[0:-6]

        centroid = alignement.centroidKey
        with open(current_file + "\\TestFiles\\GetSequenceCentroid\\" + test_case, 'r') as f:
            centroid_file = f.read()
            assert centroid == centroid_file

        aligned = alignement.aligned
        for key in aligned.keys():
            expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\AlignSequence\\" + test_case + "\\" + key, '.fasta')
            assert aligned[key] == expected

        starAlignement = alignement.heuristicMSA            
        expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\StarAlignement\\" + test_case, '.fasta')
        assert starAlignement == expected

        windowed = alignement.windowed
        for key in windowed.keys():
            expected = AlignSequences.fileToDict(current_file + "\\TestFiles\\SlidingWindow\\" + test_case + "\\" + key, '.fasta')
            assert windowed[key] == expected

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

    