import pandas as pd
import json
from aphylogeo import utils
from aphylogeo.params import Params
from Bio import Phylo
from io import StringIO
from pathlib import Path

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

# Compare trees and filter best matches
print(f"[INFO] Comparing and filtering phylogenetic trees...")
filtered_results = utils.filterResults(climatic_trees, genetic_trees, df)

# Output path for final results
output_dir = Path("results/cumacea")
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "output_cumacea.csv"

print(f"[INFO] Writing filtered results to {output_path}")
utils.writeOutputFile(filtered_results, output_path)
