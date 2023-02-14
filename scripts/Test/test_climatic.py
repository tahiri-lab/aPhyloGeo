import aPhyloGeo
import os
import pandas as pd

current_file = os.path.dirname(__file__)
climaticDataFilePath = '../datasets/5seq/geo.csv'
genetic_test_cases = ['seq small.fasta', 'seq very small.fasta', 'seq.fasta']
climatic_test_cases = ['geo.csv']

def test_openCSV():

    df1 = pd.read_csv(climaticDataFilePath)
    df2 = aPhyloGeo.openCSV(climaticDataFilePath)
    assert df1.equals(df2)


def test_climaticPipeline():
    dir = current_file + '\\TestFiles\\Dissimilarity'
    
    for test_case in climatic_test_cases:
        curent_dir = dir + '\\' + test_case
        dir_tree = current_file + '\\TestFiles\\CreateTree\\' + test_case
        df = aPhyloGeo.openCSV(climaticDataFilePath)
        trees = {}

        for filename in os.listdir(curent_dir):
            if filename.endswith(".csv"):
                filename_without_ext = filename[0:-4]

                with open(curent_dir + '\\' + filename, 'r') as matrix_expected:
                    matrix = aPhyloGeo.getDissimilaritiesMatrix(df, 'id' , filename_without_ext)    
                    # test getDissimilaritiesMatrix
                    assert str(matrix) == matrix_expected.read()

                    trees[filename_without_ext] = aPhyloGeo.createTree(matrix)

                    with open(dir_tree + '\\' + filename_without_ext + '.txt', 'r') as create_tree_expected:

                        # test createTree
                        assert str(aPhyloGeo.createTree(matrix)) == create_tree_expected.read()

        # test leastSquare
        assert aPhyloGeo.leastSquare(trees['ALLSKY_SFC_SW_DWN'],trees['T2M']) == 2.1550089999999997

