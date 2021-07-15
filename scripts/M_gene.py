from pipeline import *

if __name__ == '__main__':
    gene = 'M'
    getGene(gene, 'ATG[GT]CAGATT(.*)TGTACAGTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
