import pipeline
import os
import subprocess

def getReferenceTree(window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    gene = 'reference'
    # pipeline.changeNameSequences()
    # number_seq = pipeline.alignSequences(gene)
    # pipeline.slidingWindow(window_size, step_size)
    files = os.listdir("output/windows")
    for file in files:
        os.system("cp output/windows/" + file + " infile")
        pipeline.createBoostrap(gene)
        pipeline.createDistanceMatrix(gene)
        pipeline.createUnrootedTree(gene)
        pipeline.createConsensusTree(gene)
        bootstrap_average = pipeline.calculateAverageBootstrap()
        # on ne tient pas en consideration l'arbre cree, on le supprime
        if bootstrap_average < float(bootstrap_threshold):
            subprocess.call(["rm", "outtree"])
        else: 
            for tree in data_names:
                pipeline.calculateRfDistance(tree)
                rfn = pipeline.standardizedRfDistance(18)
                if rfn <= rf_threshold:
                    # dans ce cas, on cree les arbres par raxml
                    pipeline.runRaxML(file, gene, tree)
                else:
                    os.system("rm outfile")



