import sys

import pandas as pd
import time
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo import utils
from aphylogeo.genetic_trees import GeneticTrees
import typer
import os

# from aphylogeo.utils import climaticPipeline, geneticPipeline, filterResults, loadSequenceFile

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

    # genetic pipeline
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
        geneticTrees = utils.geneticPipeline(alignements.msa)
        trees = GeneticTrees(trees_dict=geneticTrees, format="newick")
        end_time = time.time()
        elapsed_time = round(end_time - start_time, 3)
        print(f"Elapsed time: {elapsed_time} seconds")

    # climatic sequence
    if climatic_tree is not None and os.path.exists(climatic_tree):
        climaticTrees = utils.load_climatic_trees(climatic_tree)
        climatic_data = utils.reverse_climatic_pipeline(climaticTrees)
    else:
        climatic_data = pd.read_csv(Params.file_name)
        climaticTrees = utils.climaticPipeline(climatic_data)

    filtered_results = utils.filterResults(climaticTrees, geneticTrees, climatic_data)

    utils.writeOutputFile(filtered_results, output)

    # save results
    if alignements is not None:
        alignements.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")

    trees.save_trees_to_json("./results/geneticTrees.json")

if __name__ == "__main__":
    app()
