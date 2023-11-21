import os
from pathlib import Path

import pandas as pd
import pytest

from aphylogeo import utils
from aphylogeo.params import Params

current_file = Path(os.path.dirname(__file__))
df = pd.read_csv("./datasets/example/geo.csv")
climatic_test_cases = ["geo.csv"]
dir = Path(current_file / "testFiles/dissimilarity")


@pytest.fixture(scope="module")
def climaticTreesSetup():
    """
    This fixture is used to run the climatic pipeline that creates the climatic trees.
    Returns:
        climaticTrees (dict): A dictionary containing the climatic trees.
    """
    Params.load_from_file(os.path.join(os.path.dirname(__file__), "pairwise_align.yaml"))
    return utils.climaticPipeline(pd.read_csv(Params.file_name))


def test_dissimilaritiesMatrix():
    """
    This test is used to test the dissimilarities matrix function.
    It contains the following tests:
        - getDissimilaritiesMatrix
    """

    print("Begin test_dissimilaritiesMatrix...")
    
   
    for test_case in climatic_test_cases:
        curent_dir = Path(dir / test_case)

        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]
                file = Path(curent_dir / filename)                
                actual_matrix = utils.getDissimilaritiesMatrix(df, "id", filename_without_ext)
                print(f"Actual matrix: {actual_matrix}")
                with open(file, "r") as expected_file:
                    expected_matrix = expected_file.read()
                print(f"file: {file}")
                print(f"Expected matrix: {expected_matrix}")
                # test getDissimilaritiesMatrix
                assert str(actual_matrix) == expected_matrix


def test_climaticPipeline():
    """
    This test is used to test the climaticPipeline function.
    It contains the following tests:
        - createTree
    """

    print("Begin test_climaticPipeline...")
    
   
    for test_case in climatic_test_cases:
        curent_dir = Path(dir / test_case)
        dir_test_case = Path(current_file / "testFiles/createTree" / test_case)

        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]                
                actual_matrix = utils.getDissimilaritiesMatrix(df, "id", filename_without_ext)
                print(f"Actual matrix: {actual_matrix}")
                file = Path(dir_test_case / (filename_without_ext + ".txt"))
                print(f"file: {file}")
                actual_tree = utils.createTree(actual_matrix)
                print(f"Actual tree: {actual_tree}")
                with open(file, "r") as expected_create_tree:
                    expected_tree = expected_create_tree.read()
                print(f"Expected tree: {expected_tree}")
                # test createTree
                assert str(actual_tree) == expected_tree


def test_createClimaticList(climaticTreesSetup):
    """
    This test is used to test the createClimaticList function.
    """

    print("Begin test_createGeneticList...")


    for test_case in climatic_test_cases:
        actual_list = utils.createClimaticList(climaticTreesSetup)
        file = Path(current_file / "testFiles/createClimaticList" / test_case / "list.txt")
        with open(file, "r") as f:
            expected_list = f.read()
            assert str(actual_list) == expected_list


def test_leastSquare():
    """
    This test is used to test the Least Square distance calculation function.
    It contains the following tests:
        - leastSquare
    """

    print("Begin test_leastSquare...")
    
   
    for test_case in climatic_test_cases:
        curent_dir = Path(dir / test_case)
        trees = {}
        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]                
                actual_matrix = utils.getDissimilaritiesMatrix(df, "id", filename_without_ext)
                trees[filename_without_ext] = utils.createTree(actual_matrix)

        # test leastSquare
        expected_least_square = 2.1550089999999997
        assert utils.leastSquare(trees["ALLSKY_SFC_SW_DWN"], trees["T2M"]) == expected_least_square


def test_robinsonFoulds():
    """
    This test is used to test the Robinson Foulds distance calculation function.
    It contains the following tests:
        - robinsonFoulds
    """

    print("Begin test_robinsonFoulds...")
    

    for test_case in climatic_test_cases:
        curent_dir = Path(dir / test_case)
        trees = {}
        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]                
                actual_matrix = utils.getDissimilaritiesMatrix(df, "id", filename_without_ext)
                trees[filename_without_ext] = utils.createTree(actual_matrix)

        # test robinsonFlouds
        expected_robinson_foulds = 4
        expected_robinson_foulds_Norm = 1.0
        assert utils.robinsonFoulds(trees["ALLSKY_SFC_SW_DWN"], trees["T2M"]) == (expected_robinson_foulds, expected_robinson_foulds_Norm)


def test_euclideanDist():
    """
    This test is used to test the Euclidean Bipartion distance calculation function.
    It contains the following tests:
        - euclideanDist
    """

    print("Begin test_euclideanDist...")
    

    for test_case in climatic_test_cases:
        curent_dir = Path(dir / test_case)
        trees = {}
        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]                
                actual_matrix = utils.getDissimilaritiesMatrix(df, "id", filename_without_ext)
                trees[filename_without_ext] = utils.createTree(actual_matrix)

         # test Euclidean
        expected_euclidean_dist = 0.5831890302466259
        assert utils.euclideanDist(trees["ALLSKY_SFC_SW_DWN"], trees["T2M"]) == expected_euclidean_dist        