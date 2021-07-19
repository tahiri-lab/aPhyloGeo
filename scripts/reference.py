import pipeline
import os
import subprocess

def getReferenceTree(window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    gene = 'reference'
    pipeline.changeNameSequences()
    number_seq = pipeline.alignSequences(gene)
    pipeline.slidingWindow(window_size, step_size)
    files = os.listdir("output/windows")
    for file in files:
        os.system("cp output/windows/" + file + " infile")
        pipeline.createBoostrap(gene)
        pipeline.createDistanceMatrix(gene)
        pipeline.createUnrootedTree(gene)
        pipeline.createConsensusTree(gene)
        bootstrap_average = pipeline.calculateAverageBootstrap()
        # on ne tient pas en consideration l'arbre cree, on le supprime
        if bootstrap_average >= float(bootstrap_threshold):
            subprocess.call(["rm", "outtree"])
            continue
        else: 
            pipeline.calculateRfDistance(data_names)
            rf = pipeline.standardizedRfDistance(number_seq)
            if rf >= rf_threshold:
                # dans ce cas, on cree les arbres par raxml
                pipeline.runRaxML(file, gene)
            else:
                os.system("rm outfile")
                continue



