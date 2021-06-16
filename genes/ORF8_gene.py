from gene import *
from pipeline import *

if __name__ == '__main__':
    gene = 'ORF8'
    getGene(gene, 'ATGAAATTTCTTGTTTT(.*)TTT[TC]ATCTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
