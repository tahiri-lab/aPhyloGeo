from pipeline import *

if __name__ == '__main__':
    gene = 'E'
    getGene(gene, 'ATGTACTCAT(.*)TCTGGTCTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
