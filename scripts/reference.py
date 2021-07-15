import pipeline

def getReferenceTree():
    gene = 'reference'
    pipeline.changeNameSequences()
    pipeline.alignSequences(gene)
    pipeline.createBoostrap(gene)
    pipeline.createDistanceMatrix(gene)
    pipeline.createUnrootedTree(gene)
    pipeline.createConsensusTree(gene)

if __name__ == '__main__':
    get_reference_tree()
