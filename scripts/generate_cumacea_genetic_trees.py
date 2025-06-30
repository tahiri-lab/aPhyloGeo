import os
import json
from io import StringIO
from Bio import AlignIO
from aphylogeo.params import Params
from aphylogeo.preprocessing import preprocess_windowed_alignment
from aphylogeo.utils import convert_alignment_to_simple_format, geneticPipeline
from aphylogeo.genetic_trees import GeneticTrees

Params.load_from_file("tests/cumacea_align.yaml")

# Input, output paths
input_path = "results/cumacea/aligned_sequences_cumacea.fasta.json"
filtered_path = "results/cumacea/filtered_aligned_sequences_cumacea.fasta.json"
formatted_path = "results/cumacea/filtered_formated_aligned_sequences_cumacea.fasta.json"
trees_path = "results/cumacea/geneticTrees_cumacea.json"

convert_alignment_to_simple_format(
    input_path=input_path,
    output_path= "results/cumacea/aligned_sequences_cumacea_formated.fasta.json"
)

# 1. Apply the genetic preprocessing by window
preprocess_windowed_alignment(
    input_path=input_path,
    threshold=Params.preprocessing_threshold_genetic,
    output_path=filtered_path
)
print(f"[✓] Preprocessing file saved in: {filtered_path}")

# Transform to simple format
convert_alignment_to_simple_format(
    input_path=filtered_path,
    output_path=formatted_path
)
print(f"[✓] Simple format file saved in: {formatted_path}")

# 3. Read the alignment in Biopython format
with open(filtered_path) as f:
    data = json.load(f)
    msaSet = {
        window: AlignIO.read(StringIO(fasta_str), "fasta")
        for window, fasta_str in data["msa"].items()
    }

# 4. Generating phylogenetics trees
geneticTrees = geneticPipeline(msaSet)

# 5. Saving the trees in Newick format
trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
trees.save_trees_to_json(trees_path)
print(f"[✓] Trees saved in: {trees_path}")
