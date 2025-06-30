import os
import pandas as pd
from aphylogeo import utils
from aphylogeo.params import Params
from aphylogeo.preprocessing import filter_low_variance_features
from pathlib import Path
import csv

Params.load_from_file("tests/cumacea_align.yaml")

#THIS ONE

specimen_col = "Sampleid"
test_dir = Path("tests/testFiles")
dataset_name = "Cumacea.csv"

preprocess_climate = int(Params.preprocessing_climatic)
thresh_climate = float(Params.preprocessing_threshold_climatic)
id_column = Params.specimen

# climatic preprocessing
if preprocess_climate:
    print(" Climatic preprocessing...")
    climatic_data = filter_low_variance_features(
        Params.file_name,
        threshold=thresh_climate,
        id_column=id_column
    )
    filtered_path = os.path.splitext(Params.file_name)[0] + "_filtered.csv"
    climatic_data.to_csv(filtered_path, index=False)
    print(f"Filtered CSV saved to: {filtered_path}")
else:
    climatic_data = pd.read_csv(Params.file_name)
climaticTrees = utils.climaticPipeline(climatic_data)

# 1. Dissimilarity and trees
dissim_dir = test_dir / "dissimilarity" / dataset_name
tree_dir = test_dir / "createTree" / dataset_name
dissim_dir.mkdir(parents=True, exist_ok=True)
tree_dir.mkdir(parents=True, exist_ok=True)

variables = [col for col in climatic_data.columns if col != specimen_col]
trees = {}

for var in variables:
    matrix = utils.getDissimilaritiesMatrix(climatic_data, specimen_col, var)
    tree = utils.createTree(matrix)
    trees[var] = tree
    if hasattr(matrix, "names") and hasattr(matrix, "matrix"):
        climatic_data_matrix = pd.DataFrame(matrix.matrix, index=matrix.names, columns=matrix.names)
        climatic_data_matrix.to_csv(dissim_dir / f"{var}.csv")
    else:
        matrix.to_csv(dissim_dir / f"{var}.csv")
    with open(tree_dir / f"{var}.txt", "w") as f:
        f.write(str(tree))

# 2. Variable List
climatic_list_dir = test_dir / "createClimaticList" / dataset_name
climatic_list_dir.mkdir(parents=True, exist_ok=True)
climaticTrees = utils.climaticPipeline(climatic_data)
actual_list = utils.createClimaticList(climaticTrees)
with open(climatic_list_dir / "list.txt", "w") as f:
    f.write(str(actual_list))

# 3. Least Square, Robinson-Foulds, Euclidean
ls_path = test_dir / "leastSquare" / dataset_name
rf_path = test_dir / "robinsonFoulds" / dataset_name
eu_path = test_dir / "euclideanDist" / dataset_name
ls_path.mkdir(parents=True, exist_ok=True)
rf_path.mkdir(parents=True, exist_ok=True)
eu_path.mkdir(parents=True, exist_ok=True)

with open(ls_path / "leastsquare.csv", "w", newline='') as lsfile, \
     open(rf_path / "robinsonfoulds.csv", "w", newline='') as rffile, \
     open(eu_path / "euclideandist.csv", "w", newline='') as eufile:
    ls_writer = csv.writer(lsfile)
    rf_writer = csv.writer(rffile)
    eu_writer = csv.writer(eufile)
    for t1 in trees:
        for t2 in trees:
            if t1 != t2:
                ls_writer.writerow([t1, t2, utils.leastSquare(trees[t1], trees[t2])])
                rf = utils.robinsonFoulds(trees[t1], trees[t2])
                rf_writer.writerow([t1, t2, rf[0], rf[1]])
                eu_writer.writerow([t1, t2, utils.euclideanDist(trees[t1], trees[t2])])