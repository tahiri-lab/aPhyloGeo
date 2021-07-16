import pipeline

def getReferenceTree():
    gene = 'reference'
    pipeline.changeNameSequences()
    pipeline.alignSequences(gene)
    pipeline.createBoostrap(gene)
    pipeline.createDistanceMatrix(gene)
    pipeline.createUnrootedTree(gene)
    pipeline.createConsensusTree(gene)
