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

def alignSequences(gene):
    sequences_file_name = gene + '_gene.fasta'
    directory_name = gene + '_gene'
    file_path = os.path.join('output', directory_name, sequences_file_name)
    subprocess.run(["./exec/muscle", "-in", file_path, "-phyiout", "infile", "-maxiters", "1", "-diags"])


def createBoostrap(gene):
    aligned_gene_file_name = 'aligned_' + gene + '_gene'
    directory_name = gene + '_gene'
    input_file_path = os.path.join('output', directory_name, aligned_gene_file_name)
    subprocess.run("./exec/seqboot")
    subprocess.run(["mv", "infile", input_file_path])
    subprocess.run(["mv", "outfile", "infile"])


def createDistanceMatrix(gene):
    bootstrap_file_name = 'bootstrap_' + gene + '_gene'
    directory_name = gene + '_gene'
    input_file_path = os.path.join('output', directory_name, bootstrap_file_name)
    subprocess.run("./exec/dnadist")
    subprocess.run(["cp", "infile", input_file_path])
    subprocess.run(["mv", "outfile", "infile"])

def createUnrootedTree(gene):
    distance_matrix_file_name = 'distance_matrix_' + gene + '_gene'
    directory_name = gene + '_gene'
    output_file_name = 'unrooted_tree_'+ gene + '_gene'
    input_file_path = os.path.join('output', directory_name, distance_matrix_file_name)
    output_file_path = os.path.join(
        'output', directory_name, output_file_name)
    subprocess.run("./exec/neighbor")
    subprocess.run(["mv", "infile", input_file_path])
    subprocess.run(["mv", "outtree", "intree"])
    subprocess.run(["mv", "outfile", output_file_path])


def createConsensusTree(gene):
    unrooted_tree_data_file_name = 'unrooted_tree_data_'+ gene + '_gene'
    outtree_file_name = 'consensus_tree_' + gene + '_gene'
    output_file_name = 'consensus_tree_data_' + gene + '_gene'
    directory_name = gene + '_gene'
    intree_file_path = os.path.join(
        'output', directory_name, unrooted_tree_data_file_name)
    outtree_file_path = os.path.join('output', directory_name, outtree_file_name)
    output_file_path = os.path.join('output', directory_name, output_file_name)
    subprocess.run("./exec/consense")
    subprocess.run(["mv", "intree", intree_file_path])
    subprocess.run(["mv", "outtree", output_file_path])
    subprocess.run(["mv", "outfile", outtree_file_path])
