import os
import pandas as pd
from aPhyloGeo.Params import Params
from aPhyloGeo import aPhyloGeo
from pathlib import Path
import pytest

current_file = Path(os.path.dirname(__file__))
climaticDataFilePath = '../datasets/5seq/geo.csv'
climatic_test_cases = ['geo.csv']


@pytest.fixture(scope="module")
def climaticTreesSetup():
    '''
    This fixture is used to run the climatic pipeline that creates the climatic trees.
    Returns:
        climaticTrees (dict): A dictionary containing the climatic trees.
    '''
    p = Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml"))
    return aPhyloGeo.climaticPipeline(p.file_name, p.names)


def test_openCSV():
    '''
    This test is used to test the openCSV function.
    '''
    
    print("Begin test_createGeneticList...")

    df1 = pd.read_csv(climaticDataFilePath)
    df2 = aPhyloGeo.openCSV(climaticDataFilePath)

    assert df1.equals(df2)


def test_climaticPipeline():
    '''
    This test is used to test the climaticPipeline function.
    It contains the following tests:
        - getDissimilaritiesMatrix
        - createTree
        - leastSquare
    '''

    print("Begin test_climaticPipeline...")
    
    dir = Path(current_file / 'TestFiles/Dissimilarity')
    
    for test_case in climatic_test_cases:

        curent_dir = Path(dir / test_case)
        dir_test_case = Path(current_file / 'TestFiles/CreateTree' / test_case)

        df = aPhyloGeo.openCSV(climaticDataFilePath)
        trees = {}

        for filename in os.listdir(curent_dir):

            if filename.endswith(".csv"):

                filename_without_ext = filename[0:-4]
                file = Path(curent_dir / filename)
                actual_matrix = aPhyloGeo.getDissimilaritiesMatrix(df, 'id', filename_without_ext) 

                with open(file, 'r') as expected_file:  
                    expected_matrix = expected_file.read()

                # test getDissimilaritiesMatrix
                assert str(actual_matrix) == expected_matrix

                trees[filename_without_ext] = aPhyloGeo.createTree(actual_matrix)
                file2 = Path(dir_test_case / (filename_without_ext + '.txt'))

                with open(file2, 'r') as expected_create_tree:
                    # test createTree
                    assert str(aPhyloGeo.createTree(actual_matrix)) == expected_create_tree.read()

        # test leastSquare
        expected_least_square = 2.1550089999999997
        assert aPhyloGeo.leastSquare(trees['ALLSKY_SFC_SW_DWN'], trees['T2M']) == expected_least_square


def test_createClimaticList(climaticTreesSetup):
    '''
    This test is used to test the createClimaticList function.
    '''

    print("Begin test_createGeneticList...")

    for test_case in climatic_test_cases:
        actual_list = aPhyloGeo.createClimaticList(climaticTreesSetup)
        with open(Path(current_file / 'TestFiles/CreateClimaticList' / test_case / "list.txt"), 'r') as f:
            expected_list = f.read()
            assert str(actual_list) == expected_list
