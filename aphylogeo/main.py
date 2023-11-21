import pandas as pd
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo import utils
from aphylogeo.genetic_trees import GeneticTrees

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

if __name__ == "__main__":
    print(titleCard + "\n")

    # geneticTrees = GeneticTrees.load_trees_from_file("./results/geneticTreesTest.json")
    # loaded_seq_alignment = Alignment.load_from_json("./results/aligned_sequences.json")

    Params.load_from_file()
    sequenceFile = utils.loadSequenceFile(Params.reference_gene_filepath)
    align_sequence = AlignSequences(sequenceFile)
    alignements = align_sequence.align()

    geneticTrees = utils.geneticPipeline(alignements.msa)
    trees = GeneticTrees(trees_dict=geneticTrees, format="newick")

    df = pd.read_csv(Params.file_name)
    climaticTrees = utils.climaticPipeline(df)
    utils.filterResults(climaticTrees, geneticTrees, df)
    alignements.save_to_json(f"./results/aligned_{Params.reference_gene_file}.json")
    trees.save_trees_to_json("./results/geneticTrees.json")
