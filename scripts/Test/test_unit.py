import aPhyloGeo
import os
import pandas as pd

current_file = os.path.dirname(__file__)
climaticDataFilePath = '../datasets/5seq/geo.csv'
test_cases = ['small.fasta']

# test climaticPipeline
def test_climaticPipeline():

    df1 = pd.read_csv(climaticDataFilePath)
    df2 = aPhyloGeo.openCSV(climaticDataFilePath)
    assert df1.equals(df2)



def test_getDissimilaritiesMatrix():
    path_to_test_files = current_file + '\\TestFiles\\Dissimilarity'
    print(path_to_test_files)
    for test_case in test_cases:
        df = aPhyloGeo.openCSV(climaticDataFilePath)
        file = path_to_test_files + '\\' + test_case

        for filename in os.listdir(file):
            if filename.endswith(".txt"):
                # open the .txt file and read it
                with open(file + '\\' + filename, 'r') as f:
                    print("Testing " + filename + " ...")
                    assert f.read() == str(aPhyloGeo.getDissimilaritiesMatrix(df, 'id' , filename[0:-4]))

   