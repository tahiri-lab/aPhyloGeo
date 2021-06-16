from gene import *
from pipeline import *

if __name__ == '__main__':
    gene = 'ORF7a'
    getGene(gene, 'ATGAAAATTAT(.*)GACAGAATGA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
