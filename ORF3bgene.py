from gene import *
from pipeline import *

if __name__ == '__main__':
    gene = 'ORF3b'
    getGene(gene, 'ATGAGGCTTT(.*)GCCTTTGTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
