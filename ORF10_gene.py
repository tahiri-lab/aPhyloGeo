from pipeline import *

if __name__ == '__main__':
    gene = 'ORF10'
    getGene(gene, 'ATGGGCTATA(.*)TCTCACATAG')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
