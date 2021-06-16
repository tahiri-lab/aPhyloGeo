from gene import *
from pipeline import *

if __name__ == '__main__':
    gene = 'N'
    getGene(gene, 'ATGTCT[CG][AT][TA]AAT(.*)TCAGGCCTAA')
    alignSequences(gene)
    createBoostrap(gene)
    createDistanceMatrix(gene)
    createUnrootedTree(gene)
    createConsensusTree(gene)
