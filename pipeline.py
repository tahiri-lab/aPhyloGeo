import subprocess
import os


# def fetchingSequences():
#     subprocess.run("./fetch_data.sh")

# def changeNameSequences():
#     sequences_file = open("output/sequences.fasta", "r")
#     list_of_lines = sequences_file.readlines()
#     for index in range(len(list_of_lines)):
#         if list_of_lines[index].startswith(">"):
#             splitted_line = list_of_lines[index].split("/")
#             name = ">" + splitted_line[2] + "\n"
#             list_of_lines[index] = name

#     sequences_file = open("output/sequences.fasta", "w")
#     sequences_file.writelines(list_of_lines)
#     sequences_file.close()

def alignSequences():
    subprocess.run(["./exec/muscle", "-in", "output/sequences.fasta", "-phyiout", "infile", "-maxiters", "1", "-diags"])
    subprocess.run(["cp", "infile", "output/aligned_sequences"])

def createBoostrap():
    subprocess.run("./exec/seqboot")
    subprocess.run(["cp", "outfile", "output/bootstrap_sequences"])
    subprocess.run(["mv", "outfile", "infile"])


def createDistanceMatrix():
    subprocess.run("./exec/dnadist")
    subprocess.run(["cp", "outfile", "output/distance_matrix"])
    subprocess.run(["mv", "outfile", "infile"])

def createUnrootedTree():
    subprocess.run("./exec/neighbor")
    subprocess.run(["cp", "outtree", "output/unrooted_tree_data"])
    subprocess.run(["mv", "outtree", "intree"])
    subprocess.run(["cp", "outfile", "output/unrooted_tree"])
    subprocess.run(["rm", "outfile"])


def createConsensusTree():
    subprocess.run("./exec/consense")
    subprocess.run(["mv", "outfile", "consensus_tree"])
    subprocess.run(["mv", "consensus_tree", "./output"])
    subprocess.run(["mv", "outtree", "consensus_tree_data"])
    subprocess.run(["mv", "consensus_tree_data", "./output"])

def cleanUpDirectory():
    subprocess.run(["rm", "infile", "intree"])


if __name__ == '__main__':
    # fetchingSequences()
    # changeNameSequences()
    alignSequences()
    createBoostrap()
    createDistanceMatrix()
    createUnrootedTree()
    createConsensusTree()
    cleanUpDirectory()
