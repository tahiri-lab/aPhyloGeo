import sys

import pandas as pd
import time
import json
from scipy.spatial.distance import pdist, squareform

from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo import utils
from aphylogeo.genetic_trees import GeneticTrees
from aphylogeo.preprocessing import filter_low_variance_features, preprocess_windowed_alignment
from aphylogeo.utils import convert_alignment_to_simple_format, get_patristic_distance_matrix, run_procrustes_analysis, run_protest_test
from Bio import AlignIO
from io import StringIO
import typer
import os
import warnings

# from aphylogeo.utils import climaticPipeline, geneticPipeline, filterResults, loadSequenceFile
warnings.filterwarnings("ignore", "'cgi' is deprecated", DeprecationWarning)
print(
    "Note: If you see a warning about 'cgi' being deprecated, don't worry! "
    "This is due to a dependency and does not affect the functionality of this tool."
)
titleCard = r"""
        ____    __               ___           ____
       /\  _`\ /\ \             /\_ \         /\  _`\
   __  \ \ \L\ \ \ \___   __  __\//\ \     ___\ \ \L\_\     __    ___
 /'__`\ \ \ ,__/\ \  _ `\/\ \/\ \ \ \ \   / __`\ \ \L_L   /'__`\ / __`\
/\ \L\.\_\ \ \/  \ \ \ \ \ \ \_\ \ \_\ \_/\ \L\ \ \ \/, \/\  __//\ \L\ \
\ \__/.\_\\ \_\   \ \_\ \_\/`____ \/\____\ \____/\ \____/\ \____\ \____/
 \/__/\/_/ \/_/    \/_/\/_/`/___/> \/____/\/___/  \/___/  \/____/\/___/
                              /\___/
                              \/__/
"""  # https://patorjk.com/software/taag/#p=display&f=Larry%203D&t=aphylogeo%20

app = typer.Typer(
    invoke_without_command=False,
    help="A tool for processing climatic and genetic data to generate phylogenetic trees.",
    callback=lambda: typer.echo(typer.style(titleCard, fg=typer.colors.GREEN))
)

@app.command()
def climate_pipeline(
        file_name: str = typer.Option(Params.PARAMETER_KEYS["file_name"], help="The name of the file containing the climatic data."),
        output: str = typer.Option("./datasets/example/climaticTrees.nwk", help="The name of the file to save the climatic trees."),
    ):
    """
    Run the climatic pipeline that process the climatic trees.

    Args:
        file_name (str): The name of the file containing the climatic data.
        output (str): The name of the file to save the climatic trees.
    """
    Params.load_from_file()

    preprocess_climate = int(Params.preprocessing_climatic)
    thresh_climate = float(Params.preprocessing_threshold_climatic)
    id_column = Params.specimen

    # Climatic preprocessing
    if preprocess_climate:
        print(" Climatic preprocessing...")
        climatic_data = filter_low_variance_features(
            file_name,
            threshold=thresh_climate,
            id_column=id_column
        )
        filtered_path = os.path.splitext(file_name)[0] + "_filtered.csv"
        climatic_data.to_csv(filtered_path, index=False)
        print(f"Filtered CSV saved to: {filtered_path}")
    else:
        climatic_data = pd.read_csv(file_name)

    climaticTrees = utils.climaticPipeline(climatic_data)

    try:
        utils.save_climatic_trees(climaticTrees, output)
    except Exception as e:
        print(f"Error saving the file: {e}")

@app.command()
def genetic_pipeline(
        reference_gene_filepath: str = typer.Option(os.path.join(Params.PARAMETER_KEYS["reference_gene_dir"], Params.PARAMETER_KEYS["reference_gene_file"]) , help="The path to the reference gene file."),
        output: str = typer.Option("./datasets/example/geneticTrees.json", help="The name of the file to save the genetic trees."),
    ):
    """
    Run the genetic pipeline that process the genetic trees.

    Args:
        reference_gene_filepath (str): The path to the reference gene file.
        output (str): The name of the file to save the genetic trees.
    """
    Params.load_from_file()

    sequenceFile = utils.loadSequenceFile(reference_gene_filepath)
    align_sequence = AlignSequences(sequenceFile)

    print("\nStarting alignement")
    start_time = time.time()
    alignements = align_sequence.align()
    alignements.save_to_json("results/aligned_sequences.fasta.json")

    align_sequence.makeMSA(
        align_sequence.slidingWindow(align_sequence.sequences), 
    )

    convert_alignment_to_simple_format(
        msa_json_path="results/aligned_sequences.fasta.json",
        output_path="results/aligned_sequences_formated.fasta.json"
    )

    # Genetic preprocessing
    if Params.preprocessing_genetic:
        print("Genetic preprocessing...")
        preprocess_windowed_alignment(
                input_path="results/aligned_sequences.fasta.json",
                threshold=Params.preprocessing_threshold_genetic,
                output_path="results/filtered_aligned_sequences.fasta.json"
        )

        genetic_path = "results/filtered_aligned_sequences.fasta.json"

        convert_alignment_to_simple_format(
        input_path=genetic_path,
        output_path="results/filtered_formated_aligned_sequences.fasta.json"
        )
        with open(genetic_path) as f:
                genetic_data = json.load(f)
                msaSet = {}
                for window, fasta_str in genetic_data["msa"].items():
                    alignment = AlignIO.read(StringIO(fasta_str), "fasta")
                    msaSet[window] = alignment
        
        geneticTrees = utils.geneticPipeline(msaSet)
        trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
    else:
        geneticTrees = utils.geneticPipeline(alignements.msa)
        trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
    end_time = time.time()
    elapsed_time = round(end_time - start_time, 3)
    print(f"Elapsed time: {elapsed_time} seconds")

    try:
        trees.save_trees_to_json(output)
    except Exception as e:
        print(f"Error saving the file: {e}")

@app.command()
def run(
    climatic_tree: str = typer.Option(None, help="The name of the file containing the climatic trees."),
    genetic_tree: str = typer.Option(None, help="The name of the file containing the genetic trees."),
    output: str = typer.Option("./results/output.csv", help="The name of the file to save the output."),
):
    """
    Run the pipelines and process the trees and phylogeographic analyses.

    Args:
        climatic_tree (str): The name of the file containing the climatic trees.
        genetic_tree (str): The name of the file containing the genetic trees.
        output (str): The name of the file to save the output.
    """
    # geneticTrees = GeneticTrees.load_trees_from_file("./results/geneticTreesTest.json")
    # loaded_seq_alignment = Alignment.load_from_json("./results/aligned_sequences.json")

    # load params
    Params.load_from_file()
    climatic_path = Params.file_name
    genetic_path = None

    # === Genentic pipeline ===
    alignements = None

    if genetic_tree is not None and os.path.exists(genetic_tree):
        geneticTrees = GeneticTrees.load_trees_from_file(genetic_tree)
        geneticTrees = geneticTrees.trees
        trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
    else:

        sequenceFile = utils.loadSequenceFile(Params.reference_gene_filepath)
        align_sequence = AlignSequences(sequenceFile)

        print("\nStarting alignement")
        start_time = time.time()
        alignements = align_sequence.align()
        msa_path = f"./results/aligned_{Params.reference_gene_file}.json"
        os.makedirs(os.path.dirname(msa_path), exist_ok=True)
        alignements.save_to_json(msa_path)

        convert_alignment_to_simple_format(
            input_path=msa_path,
            output_path="results/aligned_sequences_formated.fasta.json"
        )

        #genetic preprocessing
        if Params.preprocessing_genetic:
            print("\n\nGenetic preprocessing...")
            
            preprocess_windowed_alignment(
                input_path="results/aligned_sequences.fasta.json",
                threshold=Params.preprocessing_threshold_genetic,
                output_path="results/filtered_aligned_sequences.fasta.json"
            )

            genetic_path = "results/filtered_aligned_sequences.fasta.json"

            convert_alignment_to_simple_format(
            input_path=genetic_path,
            output_path="results/filtered_formated_aligned_sequences.fasta.json"
            )
            with open(genetic_path) as f:
                genetic_data = json.load(f)
                msaSet = {}
                for window, fasta_str in genetic_data["msa"].items():
                    alignment = AlignIO.read(StringIO(fasta_str), "fasta")
                    msaSet[window] = alignment
            geneticTrees = utils.geneticPipeline(msaSet)
            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
        else:
            geneticTrees = utils.geneticPipeline(alignements.msa)
            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
        end_time = time.time()
        elapsed_time = round(end_time - start_time, 3)
        print(f"Elapsed time: {elapsed_time} seconds")

    # === Climatic sequence ===
    if climatic_tree is not None and os.path.exists(climatic_tree):
        climaticTrees = utils.load_climatic_trees(climatic_tree)
        climatic_data = pd.read_csv(climatic_path)
        climatic_data = utils.reverse_climatic_pipeline(climaticTrees, climatic_data)
    else:
        # load parameters
        preprocess_climate = int(Params.preprocessing_climatic)
        thresh_climate = float(Params.preprocessing_threshold_climatic)
        id_column = Params.specimen

        # climatic preprocessing
        if preprocess_climate:
            print("\n\nClimatic preprocessing...")
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

    filtered_results = utils.filterResults(climaticTrees, geneticTrees, climatic_data)
    utils.writeOutputFile(filtered_results, output)

    # === Statistical tests ===
    try:
            # Prepare climatic distance matrix
            climatic_matrix = climatic_data.drop(columns=[Params.specimen])
            climatic_dist = squareform(pdist(climatic_matrix, metric="euclidean"))

            # Prepare genetic distance matrix (simple version from trees)
            genetic_dist = get_patristic_distance_matrix(geneticTrees)
            genetic_matrix = pd.DataFrame(genetic_dist)
    except Exception as e:
        print(f"Could not compute Statistical test: {e}")

    if Params.statistical_test == '0' or Params.statistical_test == '1':
    # === Mantel test (statistical correlation) ===
        print(f"\n\nRunning Mantel test {Params.mantel_test_method} with {Params.permutations_mantel_test} permutations... ")
        # Run Mantel test
        r, p, n = utils.run_mantel_test(genetic_dist, climatic_dist, Params.permutations_mantel_test, Params.mantel_test_method)
        print(f"Mantel test result: \n r = {r:.3f} Correlation coefficient \n p = {p:.4f} Significance level \n n = {n} Number of observations\n")
    
        # Make the row for summary
        result_row = {
            "window": "Mantel test result:",
            "r (Correlation coefficient)": r,
            "p (Significance level)": p,
            "n (Number of observations)": n
        }
        # Load the input
        if os.path.exists(output):
            df_output = pd.read_csv(output)
        else:
            df_output = pd.DataFrame()

        # Transform the row to DataFrame
        summary_df = pd.DataFrame([result_row])
        df_output = pd.concat([df_output, summary_df], ignore_index=True)

        # Guardar
        df_output.to_csv(output, index=False)

    if Params.statistical_test == '0' or Params.statistical_test == '2':
        # === Procrustes analysis ===
        print("\n\nRunning Procrustes analysis...")
        m2, genetic_transf, climatic_transf = run_procrustes_analysis(genetic_matrix, climatic_matrix)
        print(f"Procrustes M² = {round(m2, 4)} (closer to 0 = better fit)")
        
        # === PROTEST ===
        print(f"\nRunning PROTEST with {Params.permutations_protest} permutations...")
        _, protest_p = utils.run_protest_test(climatic_matrix, genetic_matrix, n_permutations=Params.permutations_protest)
        print(f"PROTEST result:\n p-value = {protest_p:.4f} (lower means more significant)")
        
        # save results
        if alignements is not None:
            os.makedirs("./results", exist_ok=True)
            alignements.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")

        trees.save_trees_to_json("./results/geneticTrees.json")

        # Make the row for summary
        result_row = {
            "window": "global",
            "Procrustes_M2": m2,
            "PROTEST_pvalue": protest_p
        }

        # Load the input
        if os.path.exists(output):
            df_output = pd.read_csv(output)
        else:
            df_output = pd.DataFrame()

        # Transform the row to DataFrame
        summary_df = pd.DataFrame([result_row])
        df_output = pd.concat([df_output, summary_df], ignore_index=True)

        # Guardar
        df_output.to_csv(output, index=False)

if __name__ == "__main__":
    app()
