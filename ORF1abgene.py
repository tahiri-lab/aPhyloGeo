from test import *
from pipeline import *

if __name__ == '__main__':
    gene = 'ORF1ab'
    getGene('ORF1ab', 'ATGGAGAGCC(.*)TAACAACTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
