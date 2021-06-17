from pipeline import *

if __name__ == '__main__':
    gene = 'S'
    getGene(gene, 'ATGTTTGTTT(.*)TTACACATAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
