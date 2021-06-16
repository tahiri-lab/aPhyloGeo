from gene import *
from pipeline import *

if __name__ == '__main__':
    gene = 'ORF1ab'
    getGene(gene, 'ATGGAGAGCC(.*)TAACAACTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
