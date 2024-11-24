import glob
import os
import re
import sys
from csv import writer as csv_writer
from io import StringIO

import dendropy
import ete3
import pandas as pd
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Applications import _Fasttree
from Bio.Phylo.Consensus import bootstrap_consensus, majority_consensus
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, _DistanceMatrix

from .multiProcessor import Multi
from .params import Params


def header():
    header = ["Gene", "Phylogeographic tree", "Name of species", "Position in ASM", "Bootstrap mean"]

    if Params.distance_method == "1":
        header.extend(["Least-Square Distance"])
    elif Params.distance_method == "2":
        header.extend(["Robinson-Foulds Distance", "Normalised RF Dist"])
    elif Params.distance_method == "3":
        header.extend(["Euclidean Distance"])
    else:
        header.extend(["Least-Square Distance"])
        header.extend(["Robinson-Foulds Distance", "Normalised RF Dist"])
        header.extend(["Euclidean Distance"])
    return header


def getDissimilaritiesMatrix(df, columnWithSpecimenName, columnToSearch):
    """
    Create a list containing the names of specimens and their corresponding minimum temperatures.

    :param df: The content of the CSV file.
    :param columnWithSpecimenName: The first column containing specimen names.
    :type columnWithSpecimenName: str
    :param columnToSearch: The column to compare with the first one.
    :type columnToSearch: str

    :return: The dissimilarities matrix.
    """
    meteoData = df[columnToSearch].tolist()
    nomVar = df[columnWithSpecimenName].tolist()
    nbrSeq = len(nomVar)
    maxValue = 0
    minValue = 0

    # First loop that allow us to calculate a matrix for each sequence
    tempTab = []

    for e in range(nbrSeq):
        # A list that will contain every distances before normalisation
        tempList = []
        for i in range(nbrSeq):
            maximum = max(float(meteoData[e]), float(meteoData[i]))
            minimum = min(float(meteoData[e]), float(meteoData[i]))
            distance = maximum - minimum
            tempList.append(float("{:.6f}".format(distance)))

        # Allow to find the maximum and minimum value for the weather value and
        # then to add the temporary list in an array
        if maxValue < max(tempList):
            maxValue = max(tempList)
        if minValue > min(tempList):
            minValue = min(tempList)
        tempTab.append(tempList)

    # Calculate normalised matrix
    tabDf = pd.DataFrame(tempTab)
    if (maxValue - minValue) != 0:
        dmDf = (tabDf - minValue) / (maxValue - minValue)
    else:
        dmDf = 0.0
    dmDf = dmDf.round(6)

    matrix = [dmDf.iloc[i, : i + 1].tolist() for i in range(len(dmDf))]
    dm = _DistanceMatrix(nomVar, matrix)
    return dm


def leastSquare(tree1, tree2):
    """
    Calculate the least square distance between two trees.

    - Trees must have the same number of leaves. 
    - Leaves must all have a twin in each tree. 
    - A tree must not have duplicate leaves.

    :param tree1: The first tree to compare.
    :type tree1: distanceTree (from Biopython)
    :param tree2: The second tree to compare.
    :type tree2: distanceTree (from Biopython)

    :return: The final distance between the two trees.
    :rtype: float
    """
    ls = 0.00
    leaves = tree1.get_terminals()

    leavesName = list(map(lambda x: x.name, leaves))
    for i in leavesName:
        leavesName.pop(0)
        for j in leavesName:
            d1 = tree1.distance(tree1.find_any(i), tree1.find_any(j))
            d2 = tree2.distance(tree2.find_any(i), tree2.find_any(j))
            ls += abs(d1 - d2)
    return ls


def robinsonFoulds(tree1, tree2):
    """
    Method that calculates the robinson foulds distance between two trees.
    
    - Trees must have the same number of leaves.
    - Leaves must all have a twin in each tree.
    - A tree must not have duplicate leaves
    
    :param tree1: The first tree to compare
    :type tree1: distanceTree object from biopython converted to Newick
    :param tree2: The second tree to compare
    :type tree2: distanceTree object from biopython converted to Newick
    
    :return: The final distance between the two
    :rtype: float
    """
    rf = 0.0
    rf_norm = 0.0
    tree1_newick = ete3.Tree(tree1.format("newick"), format=1)
    tree2_newick = ete3.Tree(tree2.format("newick"), format=1)

    rf, rf_max, common_leaves, x2, x3, x4, x5 = tree1_newick.robinson_foulds(tree2_newick, unrooted_trees=True)
    if len(common_leaves) == 0:
        rf = 0

    if rf / rf_max != 0:
        rf_norm = rf / rf_max

    return rf, rf_norm


def euclideanDist(tree1, tree2):
    """
    Method that calculate the Euclidean distance (a.k.a. Felsenstein's 2004 "branch length
    distance") between two trees based on ``edge_weight_attr``.

    Trees need to share the same |TaxonNamespace| reference.
    The bipartition bitmasks of the trees must be correct for the current tree
    structures (by calling :meth:`Tree.encode_bipartitions()` method)

    :param tree1: The first tree to compare
    :type tree1: distanceTree object from biopython converted to DendroPY format Newick
    :param tree2: The second tree to compare
    :type tree2: distanceTree object from biopython converted to DendroPY format Newick
    
    :return: The final Euclidean distance between the two
    :rtype: float
    """
    ed = 0
    tns = dendropy.TaxonNamespace()
    tree1_tc = dendropy.Tree.get(data=tree1.format("newick"), schema="newick", taxon_namespace=tns)
    tree2_tc = dendropy.Tree.get(data=tree2.format("newick"), schema="newick", taxon_namespace=tns)
    tree1_tc.encode_bipartitions()
    tree2_tc.encode_bipartitions()
    ed = dendropy.calculate.treecompare.euclidean_distance(tree1_tc, tree2_tc)

    return ed


def createTree(dm):
    """
    Create a dna tree from content coming from a fasta file.

    :param dm: The content used to create the tree.
    :return: The tree created.
    """
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def climaticPipeline(df):
    """
    Creates a dictionnary containing the climatic Trees
    
    :param df: Data to use to create the trees
    :return: The climatic tree dictionnary
    :rtype: dict
    """
    names = Params.names
    trees = {}
    for i in range(1, len(names)):
        dm = getDissimilaritiesMatrix(df, names[0], names[i])
        trees[names[i]] = createTree(dm)
    return trees

def reverse_climatic_pipeline(trees):
    """
    Converts a dictionary of climatic trees back into a DataFrame.

    :param trees: The climatic trees dictionary.
    :return: The DataFrame containing the climatic data.
    :rtype: pd.DataFrame
    """
    data = []
    for key, tree in trees.items():
        for leaf in tree.get_terminals():
            data.append([leaf.name, key, leaf.branch_length])
    df = pd.DataFrame(data, columns=[Params.names[0], 'Variable', 'Value'])
    return df

def createBoostrap(msaSet: dict, bootstrapAmount):
    """
    Create a tree structure from sequences given by a dictionnary.
    
    :param msaSet: A dictionary containing the multiple sequence alignments.
    :type msaSet: dict
    :param bootstrapAmount:
    :type bootstrapAmount: int
    
    :return: A dictionary with the trees for each sequence.
    :rtype: dict
    """
    constructor = DistanceTreeConstructor(DistanceCalculator("identity"))

    # creation of intermidiary array
    # each array is a process
    array = []
    for key in msaSet.keys():
        array.append([msaSet, constructor, key, bootstrapAmount])

    # multiprocessing
    print("Creating bootstrap variations with multiplyer of : ", bootstrapAmount)

    result = Multi(array, bootSingle).processingSmallData()

    result = sorted(result, key=lambda x: int(x[1].split("_")[0]))

    # reshaping the output into a readble dictionary
    return {i[1]: i[0] for i in result}


def bootSingle(args):
    """
    Perform a bootstrap consensus analysis on a multiple sequence alignment (MSA).

    :param args: A list of arguments, which includes:
    
        - msaSet (dict): A dictionary of multiple sequence alignments.
        - constructor (function): A function used to construct trees from the MSA.
        - key (str): The key to access a specific MSA from the msaSet.
        - bootstrapAmount (int): The number of bootstrap samples to generate.
    
    :returns: 
    
        - result (Tree): The consensus tree obtained from the bootstrap analysis.
        - key (str): The key corresponding to the MSA used for this result.
        
    :rtype: list
    """
    msaSet = args[0]
    constructor = args[1]
    key = args[2]
    bootstrapAmount = args[3]

    result = bootstrap_consensus(msaSet[key], bootstrapAmount, constructor, majority_consensus)
    return [result, key]


def calculateAverageBootstrap(tree):
    """
    Calculate the average bootstrap confidence of a tree.

    :param tree: The tree from which to calculate the average confidence.

    :return: The average bootstrap (confidence) value of the tree.
    :rtype: float
    """
    leaves = tree.get_nonterminals()
    treeConfidences = list(map(lambda x: x.confidence, leaves))
    treeConfidences = [conf for conf in treeConfidences if conf is not None]

    # Convert confidences to % if FastTree is used
    if Params.tree_type == "2":
        treeConfidences = [int(conf * 100) for conf in treeConfidences]

    totalConfidence = 0
    if treeConfidences:
        for confidences in treeConfidences:
            totalConfidence += confidences
        averageBootsrap = totalConfidence / len(treeConfidences)
    else:
        averageBootsrap = 100.0

    return averageBootsrap


def fasttreeCMD(input_fasta, boot, nt):
    """
    Create the command line executed to create a phylogenetic tree using FastTree application

    :param input_fasta: The input FASTA file containing the multiple sequence alignment.
    :type input_fasta: str
    :param output: The relative path of the output tree generated upon FastTree execution.
    :type output: str
    :param boot: The value of bootstrap to be used for tree generation (default is 100).
    :type boot: int, optional
    :param nt: Whether the sequences are nucleotides (default is True for nucleotide sequences, set to False for protein sequences).
    :type nt: bool, optional
    
    :return: The FastTree command line object.
    :rtype: FastTreeCommandline
    """
    if sys.platform == "win32":
        fasttree_exe = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\FastTree.exe"
    elif (sys.platform == "linux") | (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "darwin"):
        fasttree_exe = r"aphylogeo/bin/FastTree"
    return _Fasttree.FastTreeCommandline(fasttree_exe, input=input_fasta, nt=nt, boot=boot)


def createTmpFasta(msaset):
    """
    Create fasta files from multiple sequences alignment

    :param msaset: A dictionary where the key is the window name, and the value is the MSA object.
    :type msaset: dict(str, AlignIO)
    """
    if sys.platform == "win32":
        [SeqIO.write(alignment, os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + f"\\aphylogeo\\bin\\tmp\\{window}.fasta", "fasta") for window, alignment in msaset.items()]
    elif (sys.platform == "linux") | (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "darwin"):
        [SeqIO.write(alignment, f"aphylogeo/bin/tmp/{window}.fasta", "fasta") for window, alignment in msaset.items()]

def fasttree(msaset, boot=1000, nt=True):
    """
    Create phylogenetic trees from a set of multiple sequence alignments using the FastTree application.

    This function creates temporary FASTA files to be used as input for FastTree. Since FastTree outputs `.tree` files, 
    the function reads them back during execution and removes the output tree files afterward.

    :param msaset: A dictionary containing the multiple sequence alignments (MSA).
    :type msaset: dict(str, AlignIO)
    :param boot: The number of bootstrap replicates to generate. Default is 100.
    :type boot: int, optional
    :param nt: Whether the sequences are nucleotides (True) or amino acids (False). Default is True (nucleotide sequences).
    :type nt: bool, optional

    :return: A dictionary where the key is the window name and the value is a tree in Newick format.
    :rtype: dict(str, Tree)
    """
    createTmpFasta(msaset)
    if sys.platform == "win32":
        alignments = glob.glob(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\tmp\\*.fasta")
        windows = [re.search("tmp\\\\(.+?).fasta", fasta).group(1) for fasta in alignments]
    elif (sys.platform == "linux") | (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "darwin"):
        alignments = glob.glob("aphylogeo/bin/tmp/*.fasta")
        windows = [re.search("tmp/(.+?).fasta", fasta).group(1) for fasta in alignments]

    # Sort windows and alignments ascendent order
    sorted_windows = sorted(windows, key=lambda s: int(re.search(r"\d+", s).group()))
    sortedindex = [windows.index(win) for win in sorted_windows]
    sorted_alignments = [alignments[i] for i in sortedindex]

    # Build FastTree command lines
    cmds = [fasttreeCMD(fasta, boot, nt) for fasta in sorted_alignments]

    # Execute cmd lines and create trees
    trees = {sorted_windows[i]: Phylo.read(StringIO(cmd()[0]), "newick") for i, cmd in enumerate(cmds)}
    if sys.platform == "win32":
        [os.remove(file) for file in glob.glob(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "\\aphylogeo\\bin\\tmp\\*.fasta")]  # Remove temp fasta files
    elif (sys.platform == "linux") | (sys.platform == "linux1") | (sys.platform == "linux2") | (sys.platform == "darwin"):
        [os.remove(file) for file in glob.glob("aphylogeo/bin/tmp/*.fasta")]  # Remove temp fasta files
    return trees


def createGeneticList(geneticTrees, bootstrap_threshold):
    """
    Create a list of Trees if the bootstrap Average is higher than
    the threshold

    :param geneticTrees: A dictionary containing the genetic trees.
    :type geneticTrees: dict
    :param bootstrap_threshold: 
    
    :return: 
     - geneticTrees: A sorted list of genetic trees.
     - bootstrapList: A list of the bootstrap average of each tree.
    """
    bootstrapList = []
    geneticList = []

    for key in geneticTrees:
        bootstrap_average = calculateAverageBootstrap(geneticTrees[key])
        if bootstrap_average >= bootstrap_threshold:
            bootstrapList.append(bootstrap_average)
            geneticList.append(key)
    return geneticList, bootstrapList


def createClimaticList(climaticTrees):
    """
    Create a list of climaticTrees
    
    :param climaticTrees: A dictionary containing the climatic trees.
    :type climaticTrees: dict
    
    :return: A list of the climatic trees.
    """
    climaticList = []
    for key in climaticTrees:
        climaticList.append(key)
    return climaticList


def getData(leavesName, dist, index, climaticList, bootstrap, genetic, csv_data, reference_gene_filename, dist_norm, dist2, dist3):
    """
    Get data from a csv file a various parameters to store into a list

    :param leavesName: The list of the actual leaves.
    :type leavesName: list
    :param dist: The least square distance between two trees.
    :type dist: float
    :param index: The index of the climatic tree.
    :type index: int
    :param climaticList: The list of climatic trees.
    :type climaticList: list
    :param bootstrap: The bootstrap values.
    :param genetic: The genetic tree.
    :param csv_data: The pf containing the data.
    :param reference_gene_filename: The name of the reference gene.
    :param dist_norm:
    :param dist2: 
    :param dist3:
    
    :return: A list containing the data.
    :rtype: list
    """

    for leave in leavesName:
        for i, row in csv_data.iterrows():
            if row.iloc[0] == leave:
                if Params.distance_method == "2":
                    return [reference_gene_filename, climaticList[index], leave, genetic, str(bootstrap), str(round(dist, 2)), str(dist_norm)]
                elif Params.distance_method == "0":
                    return [
                        reference_gene_filename,
                        climaticList[index],
                        leave,
                        genetic,
                        str(bootstrap),
                        str(round(dist, 2)),
                        str(dist2),
                        str(dist_norm),
                        str(round(dist3, 2)),
                    ]
                else:
                    return [reference_gene_filename, climaticList[index], leave, genetic, str(bootstrap), str(round(dist, 2))]


def writeOutputFile(data):
    """
    Write the datas from data list into a new csv file

    :param data: The list containing the final data.
    :type data: list
    """
    print("Writing the output file")
    directory = os.path.abspath("./results")
    os.makedirs(directory, exist_ok=True)
    with open("./results/output.csv", "w", encoding="UTF8") as f:
        writer = csv_writer(f)
        writer.writerow(header())
        for i in range(len(data)):
            writer.writerow(data[i])
        f.close


def filterResults(
    climaticTrees,
    geneticTrees,
    csv_data,
    create_file=True,
):
    """
    Create the final datas from the Climatic Tree and the Genetic Tree

    :param climaticTrees: The dictionary containing the climatic trees.
    :type climaticTrees: dict
    :param geneticTrees: The dictionary containing the genetic trees.
    :type geneticTrees: dict
    :param csv_data: The dataframe containing the data from the CSV file.
    :type csv_data: pd.DataFrame
    :param create_file: Whether to create a file or not. Default is True.
    :type create_file: bool
    
    :return: The final data in a dictionary format.
    :rtype: dict
    """

    # Create a list of the tree if the bootstrap is superior to the
    # bootstrap treshold
    geneticList, bootstrapList = createGeneticList(geneticTrees, Params.bootstrap_threshold)

    # Create a list with the climatic trees name
    climaticList = createClimaticList(climaticTrees)

    data = []
    # Compare every genetic trees with every climatic trees. If the returned
    # value is inferior or equal to the distance threshold,
    # we keep the data
    while len(geneticList) > 0:
        current_genetic = geneticList.pop(0)
        current_bootstrap = bootstrapList.pop(0)

        leaves = geneticTrees[current_genetic].get_terminals()
        leavesName = list(map(lambda x: x.name, leaves))

        for i in range(len(climaticTrees.keys())):
            if Params.distance_method == "1":
                ls = leastSquare(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                if ls is None:
                    raise Exception("The LS distance is not calculable" + " for {aligned_file}.")
                if ls <= Params.dist_threshold:
                    data.append(
                        getData(
                            leavesName,
                            ls,
                            i,
                            climaticList,
                            current_bootstrap,
                            current_genetic,
                            csv_data,
                            Params.reference_gene_filepath,
                            None,
                            None,
                            None,
                        )
                    )
            elif Params.distance_method == "2":
                rf, rf_norm = robinsonFoulds(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                if rf is None:
                    raise Exception("The RF distance is not calculable" + " for {aligned_file}.")
                if rf <= Params.dist_threshold:
                    data.append(
                        getData(
                            leavesName,
                            rf,
                            i,
                            climaticList,
                            current_bootstrap,
                            current_genetic,
                            csv_data,
                            Params.reference_gene_filepath,
                            rf_norm,
                            None,
                            None,
                        )
                    )
            elif Params.distance_method == "3":
                ed = euclideanDist(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                if ed is None:
                    raise Exception("The Euclidean (DendroPY) distance is not calculable" + " for {aligned_file}.")
                if ed <= Params.dist_threshold:
                    data.append(
                        getData(
                            leavesName,
                            ed,
                            i,
                            climaticList,
                            current_bootstrap,
                            current_genetic,
                            csv_data,
                            Params.reference_gene_filepath,
                            None,
                            None,
                            None,
                        )
                    )
            elif Params.distance_method == "0":
                ls = leastSquare(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                rf, rf_norm = robinsonFoulds(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                ed = euclideanDist(geneticTrees[current_genetic], climaticTrees[climaticList[i]])
                if ls is None or rf is None or ed is None:
                    raise Exception("The LS distance is not calculable" + " for {aligned_file}.")
                if ls <= Params.dist_threshold:
                    data.append(
                        getData(
                            leavesName,
                            ls,
                            i,
                            climaticList,
                            current_bootstrap,
                            current_genetic,
                            csv_data,
                            Params.reference_gene_filepath,
                            rf_norm,
                            rf,
                            ed,
                        )
                    )
            else:
                raise ValueError("Invalid distance method")

    if create_file:
        # We write the datas into an output csv file
        writeOutputFile(data)
    return format_to_csv(data)


def format_to_csv(data):
    """
    Format the data to a csv file
    
    :param data: The data to format.
    :type data: list
    
    :return: The formatted data in a dictionary format.
    :rtype: dict
    """
    result = {}

    for h in header():
        result[h] = []

    for row in data:
        for col_index in range(len(row)):
            result[header()[col_index]].append(row[col_index])

    return result


def geneticPipeline(seq_alignement):
    """
    Get the genetic Trees from the initial file datas so we
    can compare every valid tree with the climatic ones. In the
    end it calls a method that create a final csv file with all
    the data that we need for the comparison

    :param seq_alignement: The multiple sequence alignment.
    :type seq_alignement: dict
    
    :return: The genetic trees.
    :rtype: dict
    """

    # GM no more needed
    # JUST TO MAKE THE DEBUG FILES
    # if os.path.exists("./debug"):
    #     shutil.rmtree("./debug")
    # if Params.makeDebugFiles:
    #     os.mkdir("./debug")
    # JUST TO MAKE THE DEBUG FILES

    if Params.tree_type == "1":
        geneticTrees = createBoostrap(seq_alignement, Params.bootstrap_amount)
    elif Params.tree_type == "2":
        geneticTrees = fasttree(seq_alignement, Params.bootstrap_amount, True)

    return geneticTrees


def loadSequenceFile(file):
    """
    Reads the .fasta file. Extract sequence ID and sequences.
    
    :param file: The file name of a .fasta file.
    :type file: str
    
    :return: The sequences in a dictionary format.
    :rtype: dict
    """
    sequences = {}
    with open(file) as sequencesFile:
        for sequence in SeqIO.parse(sequencesFile, "fasta"):
            sequences[sequence.id] = sequence.seq
    return sequences

def save_climatic_trees(climatic_trees, file_path):
    """
    Save the climatic trees to a file.

    :param dict climatic_trees: The climatic trees to save.
    :param str file_path: The path to the file where the climatic trees will be saved.
    """

    with open(file_path, 'w') as file:
        for key, tree in climatic_trees.items():
            file.write(f">{key}\n")
            Phylo.write(tree, file, "newick")


def load_climatic_trees(file_path):
    """
    Load climatic trees from a file.

    :param str file_path: Path to the target file where the climatic trees will be loaded
    :returns: The loaded climatic trees
    :rtype: dict
    """
    climatic_trees = {}

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):
            key = lines[i].strip()[1:]
            tree = Phylo.read(StringIO(lines[i + 1]), "newick")
            climatic_trees[key] = tree
    return climatic_trees
