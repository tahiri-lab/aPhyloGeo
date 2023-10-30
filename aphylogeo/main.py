import pandas as pd
from aphylogeo.alignement import AlignSequences, Alignment

from aphylogeo.params import Params
from aphylogeo import utils

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

    sequenceFile = utils.loadSequenceFile(Params().reference_gene_file)
    seq_alignment = AlignSequences(sequenceFile).align()
    # Phylo.write(tree1, "data/tree1.nwk", "newick")
    seq_alignment.save_to_json("./debug/sequences_aligned.json")

    loaded_seq_alignment = Alignment.load_from_json("./debug/sequences_aligned.json")

    geneticTrees = utils.geneticPipeline(seq_alignment.msa)

    # Todo get trees in geneticTrees dict
    # Phylo.write(geneticTrees, "./tree1.nwk", "newick")
    df = pd.read_csv(Params().file_name)
    climaticTrees = utils.climaticPipeline(df)
    utils.filterResults(climaticTrees, geneticTrees, df)
