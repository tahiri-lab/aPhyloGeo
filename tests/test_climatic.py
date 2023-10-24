import os
import pandas as pd
from aphylogeo.params import Params
from aphylogeo import utils
from pathlib import Path
import pytest

current_file = Path(os.path.dirname(__file__))
df = pd.read_csv('./datasets/5seq/geo.csv')
climatic_test_cases = ['geo.csv']


@pytest.fixture(scope="module")
def climaticTreesSetup():
    '''
    This fixture is used to run the climatic pipeline that creates the climatic trees.
    Returns:
        climaticTrees (dict): A dictionary containing the climatic trees.
    '''
    p = Params(os.path.join(os.path.dirname(__file__), "params_very_small.yaml"))
    return utils.climaticPipeline(pd.read_csv(p.file_name), p.names)


def test_climaticPipeline():
    '''
    This test is used to test the climaticPipeline function.
    It contains the following tests:
        - getDissimilaritiesMatrix
        - createTree
        - leastSquare
    '''

    print("Begin test_climaticPipeline...")
    
    dir = Path(current_file / 'testFiles/dissimilarity')
    
    for test_case in climatic_test_cases:

        curent_dir = Path(dir / test_case)
        dir_test_case = Path(current_file / 'testFiles/createTree' / test_case)

        trees = {}

        for filename in os.listdir(curent_dir):

            if filename.endswith(".csv"):

                filename_without_ext = filename[0:-4]
                file = Path(curent_dir / filename)
                actual_matrix = utils.getDissimilaritiesMatrix(df, 'id', filename_without_ext) 
                print(f'Actual matrix: {actual_matrix}')
                with open(file, 'r') as expected_file:  
                    expected_matrix = expected_file.read()
                print(f'file: {file}')
                print(f'Expected matrix: {expected_matrix}')
                # test getDissimilaritiesMatrix
                assert str(actual_matrix) == expected_matrix

                trees[filename_without_ext] = utils.createTree(actual_matrix)
                file2 = Path(dir_test_case / (filename_without_ext + '.txt'))

                with open(file2, 'r') as expected_create_tree:
                    # test createTree
                    assert str(utils.createTree(actual_matrix)) == expected_create_tree.read()        

        # test leastSquare
        expected_least_square = 2.1550089999999997
        assert utils.leastSquare(trees['ALLSKY_SFC_SW_DWN'], trees['T2M']) == expected_least_square
       
        # test robinsonFlouds
        expected_robinson_foulds = 4
        expected_robinson_fouldsMAX = 1.0
        assert utils.robinsonFoulds(trees['ALLSKY_SFC_SW_DWN'], trees['T2M']) == (expected_robinson_foulds, expected_robinson_fouldsMAX)


def test_createClimaticList(climaticTreesSetup):
    '''
    This test is used to test the createClimaticList function.
    '''

    print("Begin test_createGeneticList...")

    for test_case in climatic_test_cases:
        actual_list = utils.createClimaticList(climaticTreesSetup)
        with open(Path(current_file / 'testFiles/createClimaticList' / test_case / "list.txt"), 'r') as f:
            expected_list = f.read()
            assert str(actual_list) == expected_list
