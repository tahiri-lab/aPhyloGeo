from pipeline import *

if __name__ == '__main__':
    gene = 'ORF6'
    getGene(gene, 'ATGTTTCATC(.*)GATTGA[CT]TAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
