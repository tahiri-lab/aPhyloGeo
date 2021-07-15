import pipeline

def getEGene():
    gene = 'E'
    pipeline.getGene(gene, 'ATGTACTCAT(.*)TCTGGTCTAA')
    pipeline.alignSequences(gene)
    pipeline.createBoostrap(gene)
    pipeline.createDistanceMatrix(gene)
    pipeline.createUnrootedTree(gene)
    pipeline.createConsensusTree(gene)
