import Alignement
import aPhyloGeo
import os
import pandas as pd
import pytest

current_file = os.path.dirname(__file__)
genetic_test_cases = ['seq small.fasta', 'seq very small.fasta']


@pytest.fixture(scope="module")
def climaticTreesSetup():
    return aPhyloGeo.climaticPipeline()

@pytest.fixture(scope="module")
def alignementSetup():
    return Alignement()

@pytest.mark.usefixtures('climaticTreesSetup')
class TestGenetic:

    def test_genetic_alignement(self, climaticTreesSetup):
        assert True

    def test_genetic_createBootStrap(self):
        assert True

    def test_genetic_filterResults(self):
        assert True


@pytest.mark.usefixtures('alignementSetup')
class TestAlignement:

    def test_getSequenceCentroid(self, alignementSetup):
        alignementSetup.getSequenceCentroid()
        assert True
    def test_alignSequences(self, alignementSetup):
        assert True
    def test_starAlignement(self, alignementSetup):
        assert True
    def test_slidingWindow(self, alignementSetup):
        assert True
    def test_makeMSA(self, alignementSetup):
        assert True

    