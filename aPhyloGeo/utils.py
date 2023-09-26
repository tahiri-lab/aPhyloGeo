import pandas as pd
import os
import shutil
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.Consensus import *
from Bio import SeqIO
from csv import writer as csv_writer
import random

from .multiProcessor import Multi
from .alignement import AlignSequences
from .params import Params

HEADER = ['Gene', 'Phylogeographic tree', 'Name of species',
          'Position in ASM', 'Bootstrap mean', 'Least-Square distance']


def getDissimilaritiesMatrix(df, columnWithSpecimenName, columnToSearch):
    """
    Creation of a list containing the names of specimens and minimums
    tempratures

    Args:
        df (content of CSV file)
        columnWithSpecimenName (first column of names)
        columnToSearch (column to compare with the first one)

    Return:
        The dissimilarities matrix

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
    dmDf = (tabDf - minValue) / (maxValue - minValue)
    dmDf = dmDf.round(6)

    matrix = [dmDf.iloc[i, :i + 1].tolist() for i in range(len(dmDf))]
    dm = _DistanceMatrix(nomVar, matrix)
    return dm


def leastSquare(tree1, tree2):
    """
    Method that calculates the least square distance between two trees.
    Trees must have the same number of leaves.
    Leaves must all have a twin in each tree.
    A tree must not have duplicate leaves
     x   x
    ╓╫╖ ╓╫╖
    123 312

    Args:
        tree1 (distanceTree object from biopython)
        tree2 (distanceTree object from biopython)

    Return:
        return result (double) the final distance between the two

    """
    ls = 0.00
    leaves = tree1.get_terminals()

    leavesName = list(map(lambda x: x.name, leaves))

    for i in leavesName:
        leavesName.pop(0)
        for j in leavesName:
            d1 = (tree1.distance(tree1.find_any(i), tree1.find_any(j)))
            d2 = (tree2.distance(tree2.find_any(i), tree2.find_any(j)))
            ls += (abs(d1 - d2))
    return ls


def drawTreesmake(trees, p):
    """
    Function that will draw the trees for each climatic variable.
    The DistanceTreeConstructor object is transformed to Newick format and
    loaded as a toytree MulTitree object. Some stylings are applied and the
    resulting trees are drawed into a .pdf in the viz/ dir.

    Args:
        trees (Dictionnary of DistanceTreeConstructor object with climatic
        variable for keys)
        p (Params object)

    """
    treesNewick = {}
    toytrees = []

    for k, v in trees.items():
        treesNewick[k] = v.format('newick')
        ttree = toytree.tree(treesNewick[k], tree_format=1)
        toytrees.append(ttree)
    mtree = toytree.mtree(toytrees)

    # Setting up the stylings for nodes
    for tree in mtree.treelist:
        tree.style.edge_align_style = {'stroke': 'black', 'stroke-width': 1}
        for node in tree.treenode.traverse():
            if node.is_leaf():
                node.add_feature('color', toytree.colors[7])
            else:
                node.add_feature('color', toytree.colors[1])
    colors = tree.get_node_values('color', show_root=1, show_tips=1)

    # Draw the climatic trees
    canvas, axes, mark = mtree.draw(nrows=round(len(mtree) / 5),
                                    ncols=len(mtree), height=400, width=1000,
                                    node_sizes=8, node_colors=colors,
                                    tip_labels_align=True)

    for i in range(len(mtree)):
        randColor = "#%03x" % random.randint(0, 0xFFF)
        axes[i].text(0, mtree.ntips, p.names[i + 1], style={'fill': randColor,
                                                            'font-size': '10px',
                                                            'font-weight': 'bold'
                                                            })

    toyplot.pdf.render(canvas, '../viz/climactic_trees.pdf')


def createTree(dm):
    '''
    Create a dna tree from content coming from a fasta file.

    Args:
        dm (content used to create the tree)

    Return:
        tree (the new tree)
    '''
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def climaticPipeline(df, names):
    '''
    Creates a dictionnary with the climatic Trees
    Args:
        df (contains the data to use)
        names (list of names of the columns to use)

    Return:
        trees (the climatic tree dictionnary)
    '''
    trees = {}
    for i in range(1, len(names)):
        dm = getDissimilaritiesMatrix(df, names[0], names[i])
        trees[names[i]] = createTree(dm)
    return trees


def createBoostrap(msaSet: dict, bootstrapAmount):
    '''
    Create a tree structure from sequences given by a dictionnary.
    Args:
        msaSet (dictionnary with multiple sequences alignment to transform into
                trees)
        bootstrapAmount
    Return:
        A dictionary with the trees for each sequence
    '''
    constructor = DistanceTreeConstructor(DistanceCalculator('identity'))

    # creation of intermidiary array
    # each array is a process
    array = []
    for key in msaSet.keys():
        array.append([msaSet, constructor, key, bootstrapAmount])

    # multiprocessing
    print("Creating bootstrap variations with multiplyer of : ",
          bootstrapAmount)

    result = Multi(array, bootSingle).processingSmallData()

    result = sorted(result, key=lambda x: int(x[1].split('_')[0]))

    # reshaping the output into a readble dictionary
    return {i[1]: i[0] for i in result}


def bootSingle(args):
    '''
    Args:
        args (list of arguments)
    Return:
        *********TO WRITE**********
    '''
    msaSet = args[0]
    constructor = args[1]
    key = args[2]
    bootstrapAmount = args[3]

    result = bootstrap_consensus(msaSet[key], bootstrapAmount, constructor,
                                 majority_consensus)
    return [result, key]


def calculateAverageBootstrap(tree):
    '''
    Calculate if the average confidence of a tree

    Args:
        tree (The tree to get the average confidence from)
    Return :
        averageBootstrap (the average Bootstrap (confidence))
    '''
    leaves = tree.get_nonterminals()
    treeConfidences = list(map(lambda x: x.confidence, leaves))
    treeConfidences.pop(0)
    totalConfidence = 0
    for confidences in treeConfidences:
        totalConfidence += confidences
    averageBootsrap = totalConfidence / len(treeConfidences)
    return averageBootsrap


def createGeneticList(geneticTrees, bootstrap_threshold):
    '''
    Create a list of Trees if the bootstrap Average is higher than
    the threshold

    Args :
        geneticTrees (a dictionnary of genetic trees)
        bootstrap_threshold
    Return :
        geneticList (a sorted list with the geneticTrees)
        bootstrapList (a list with the bootstrap average of each tree)
    '''
    bootstrapList = []
    geneticList = []

    for key in geneticTrees:
        bootstrap_average = calculateAverageBootstrap(geneticTrees[key])
        if bootstrap_average >= bootstrap_threshold:
            bootstrapList.append(bootstrap_average)
            geneticList.append(key)
    return geneticList, bootstrapList


def createClimaticList(climaticTrees):
    '''
    Create a list of climaticTrees

    Args :
        climaticTrees (a dictionnary of climatic trees)
    Return :
        climaticList(a list with the climaticTrees)
    '''
    climaticList = []
    for key in climaticTrees:
        climaticList.append(key)
    return climaticList


def getData(leavesName, ls, index, climaticList, bootstrap, genetic, csv_data, reference_gene_filename):
    '''
    Get data from a csv file a various parameters to store into a list

    Args :
        leavesName (the list of the actual leaves)
        ls (least square distance between two trees)
        climaticList (the list of climatic trees)
        bootstrap (bootstrap values)
        genetic  (genetic tree)
        csv_data (the pf containing the data)
        reference_gene_filename
    '''

    for leave in leavesName:
        for i, row in csv_data.iterrows():
            if row[0] == leave:
                return [reference_gene_filename, climaticList[index],
                        leave, genetic,
                        str(bootstrap), str(round(ls, 2))]


def writeOutputFile(data):
    '''
    Write the datas from data list into a new csv file

    Args :
        data (the list contaning the final data)
    '''
    print("Writing the output file")

    with open("output.csv", "w", encoding="UTF8") as f:
        writer = csv_writer(f)
        writer.writerow(HEADER)
        for i in range(len(data)):
            writer.writerow(data[i])
        f.close


def filterResults(climaticTrees, geneticTrees, bootstrap_threshold, ls_threshold, csv_data, reference_gene_filename, create_file=True):
    '''
    Create the final datas from the Climatic Tree and the Genetic Tree

    Args :
        climaticTrees (the dictionnary containing every climaticTrees)
        geneticTrees (the dictionnary containing every geneticTrees)
        bootstrap_threshold (the bootstrap threshold)
        ls_threshold (the least square threshold)
        csv_data (dataframe containing the data from the csv file)
        reference_gene_filename (the name of the reference gene)
    '''
    # Create a list of the tree if the bootstrap is superior to the
    # bootstrap treshold
    geneticList, bootstrapList = createGeneticList(geneticTrees, bootstrap_threshold)

    # Create a list with the climatic trees name
    climaticList = createClimaticList(climaticTrees)

    data = []
    # Compare every genetic trees with every climatic trees. If the returned
    # value is inferior or equal to the (Least-Square distance) LS threshold,
    # we keep the data
    while (len(geneticList) > 0):

        current_genetic = geneticList.pop(0)
        current_bootstrap = bootstrapList.pop(0)

        leaves = geneticTrees[current_genetic].get_terminals()
        leavesName = list(map(lambda x: x.name, leaves))

        for i in range(len(climaticTrees.keys())):
            ls = leastSquare(geneticTrees[current_genetic],
                             climaticTrees[climaticList[i]])
            if ls is None:
                raise Exception('The LS distance is not calculable' + 'pour {aligned_file}.')
            if ls <= ls_threshold:
                data.append(getData(leavesName, ls, i, climaticList,
                                    current_bootstrap, current_genetic, csv_data, reference_gene_filename))

    if create_file:
        # We write the datas into an output csv file
        writeOutputFile(data)
    return format_to_csv(data)


def format_to_csv(data):
    '''
    Format the data to a csv file

    Args:
        data: array of arrays of data

    Returns:
        result: dict with key header and value array of data
    '''
    result = {}

    for h in HEADER:
        result[h] = []

    for row in data:
        for col_index in range(len(row)):
            result[HEADER[col_index]].append(row[col_index])

    return result


def geneticPipeline(climaticTrees, csv_data, p=Params(), alignementObject=None):
    '''
    Get the genetic Trees from the initial file datas so we
    can compare every valid tree with the climatic ones. In the
    end it calls a method that create a final csv file with all
    the data that we need for the comparison

    Args:
        climaticTrees (the dictionnary of climaticTrees)
        p (the Params object)[optionnal]
        aligneementObject (the AlignSequences object)[optionnal]
    '''
    # JUST TO MAKE THE DEBUG FILES
    if os.path.exists("./debug"):
        shutil.rmtree("./debug")
    if p.makeDebugFiles:
        os.mkdir("./debug")
    # JUST TO MAKE THE DEBUG FILES

    if alignementObject is None:
        sequences = openFastaFile(p.reference_gene_file)
        alignementObject = AlignSequences(sequences, p.window_size, p.step_size, p.makeDebugFiles, p.bootstrapAmount, p.alignment_method, p.reference_gene_file)

    msaSet = alignementObject.msaSet
    geneticTrees = createBoostrap(msaSet, p.bootstrapAmount)
    return filterResults(climaticTrees, geneticTrees, p.bootstrap_threshold, p.ls_threshold, csv_data, p.reference_gene_filename)


def openFastaFile(file):
    '''
    Reads the .fasta file. Extract sequence ID and sequences.

    Args:
        file (String) the file name of a .fasta file

    Return:
        sequences (dictionnary)
    '''
    sequences = {}
    with open(file) as sequencesFile:
        for sequence in SeqIO.parse(sequencesFile, "fasta"):
            sequences[sequence.id] = sequence.seq
    return sequences
