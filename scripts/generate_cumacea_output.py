import pandas as pd
import json
from aphylogeo import utils
from aphylogeo.params import Params
from Bio import Phylo
from io import StringIO
from pathlib import Path
from aphylogeo.utils import get_patristic_distance_matrix
from scipy.spatial.distance import pdist, squareform

# Load parameters
Params.load_from_file("tests/cumacea_align.yaml")

# Load climatic data
df = pd.read_csv(Params.file_name)

# Path to genetic trees (JSON format, window -> Newick string)
genetic_tree_path = Path("results/cumacea/geneticTrees_cumacea.json")
if not genetic_tree_path.exists():
    raise FileNotFoundError(f"[ERROR] File not found: {genetic_tree_path}")

print(f"[INFO] Loading genetic trees from: {genetic_tree_path}")
with open(genetic_tree_path) as f:
    genetic_trees = json.load(f)

# Parse Newick strings to Bio.Phylo tree objects
genetic_trees = {
    k: Phylo.read(StringIO(tree_str), "newick")
    for k, tree_str in genetic_trees.items()
}

# Generate climatic trees from the dataframe using your defined pipeline
print(f"[INFO] Generating climatic trees from climatic data...")
climatic_trees = utils.climaticPipeline(df)

# === Mantel test (statistical correlation) ===
try:
    print(f"\n\nRunning Mantel test {Params.mantel_test_method} with {Params.permutations_mantel_test} permutations... ")

    # Prepare climatic distance matrix
    climatic_matrix = df.drop(columns=[Params.specimen])
    climatic_dist = squareform(pdist(climatic_matrix, metric="euclidean"))

    # Prepare genetic distance matrix (simple version from trees)
    genetic_dist = get_patristic_distance_matrix(genetic_trees)

    # Run Mantel test
    r, p, n = utils.run_mantel_test(genetic_dist, climatic_dist, Params.permutations_mantel_test, Params.mantel_test_method)

    print(f"Mantel test result: \n r = {r:.3f} Correlation coefficient \n p = {p:.4f} Significance level \n n = {n} Number of observations\n")
except Exception as e:
    print(f"Could not compute Mantel test: {e}")


# Compare trees and filter best matches
print(f"[INFO] Comparing and filtering phylogenetic trees...")
filtered_results = utils.filterResults(climatic_trees, genetic_trees, df)

# Output path for final results
output_dir = Path("results/cumacea")
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "output_cumacea.csv"

print(f"[INFO] Writing filtered results to {output_path}")
utils.writeOutputFile(filtered_results, output_path)
