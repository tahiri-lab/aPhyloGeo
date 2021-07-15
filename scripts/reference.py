from pipeline import *

def get_reference_tree():
    gene = 'reference'
    changeNameSequences()
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)

    