from pipeline import *

if __name__ == '__main__':
    gene = 'ORF3a'
    getGene(gene, 'ATGGATTTGT(.*)GCCTTTGTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
