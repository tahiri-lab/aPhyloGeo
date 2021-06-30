import subprocess
import os
import re

def changeNameSequences():
    sequences_file = open("output/reference_gene/reference_gene.fasta", "r")
    list_of_lines = sequences_file.readlines()
    for index in range(len(list_of_lines)):
        if list_of_lines[index].startswith(">"):
            splitted_line = list_of_lines[index].split("/")
            name = ">" + splitted_line[2] + "\n"
            list_of_lines[index] = name

    sequences_file = open("output/reference_gene/reference_gene.fasta", "w")
    sequences_file.writelines(list_of_lines)
    sequences_file.close()


def getGene(gene, pattern): 
    sequences_file = open("output/reference_gene/reference_gene.fasta", "r").read()
    list_of_sequences = sequences_file.split(">")
    s = pattern
    directory_name = gene + "_gene"
    file_name = gene + "_gene.fasta"
    path = os.path.join("output", directory_name, file_name)
    new_file = open(path, "w")
    for index in range(len(list_of_sequences)):
        if list_of_sequences[index] == "":
            continue
        name = list_of_sequences[index].split("\n")[0]
        gene_sequence = list_of_sequences[index].replace("\n", "")
        gene_sequence = (re.search(s, gene_sequence).group())
        new_file.writelines(">" + name + "\n")
        new_file.writelines(gene_sequence + "\n")

    new_file.close()

def alignSequences(gene):
    sequences_file_name = gene + '_gene.fasta'
    directory_name = gene + '_gene'
    file_path = os.path.join('output', directory_name, sequences_file_name)
    subprocess.call(["./exec/muscle", "-in", file_path, "-physout", "infile", "-maxiters", "1", "-diags"])


def createBoostrap(gene):
    aligned_gene_file_name = 'aligned_' + gene + '_gene'
    directory_name = gene + '_gene'
    input_file_path = os.path.join('output', directory_name, aligned_gene_file_name)
    subprocess.call("./exec/seqboot")
    subprocess.call(["mv", "infile", input_file_path])
    subprocess.call(["mv", "outfile", "infile"])


def createDistanceMatrix(gene):
    bootstrap_file_name = 'bootstrap_' + gene + '_gene'
    directory_name = gene + '_gene'
    input_file_path = os.path.join('output', directory_name, bootstrap_file_name)
    subprocess.call("./exec/dnadist")
    subprocess.call(["cp", "infile", input_file_path])
    subprocess.call(["mv", "outfile", "infile"])

def createUnrootedTree(gene):
    distance_matrix_file_name = 'distance_matrix_' + gene + '_gene'
    directory_name = gene + '_gene'
    output_file_name = 'unrooted_tree_'+ gene + '_gene'
    input_file_path = os.path.join('output', directory_name, distance_matrix_file_name)
    output_file_path = os.path.join(
        'output', directory_name, output_file_name)
    subprocess.call("./exec/neighbor")
    subprocess.call(["mv", "infile", input_file_path])
    subprocess.call(["mv", "outtree", "intree"])
    subprocess.call(["mv", "outfile", output_file_path])


def createConsensusTree(gene):
    unrooted_tree_data_file_name = 'unrooted_tree_data_'+ gene + '_gene'
    outtree_file_name = 'consensus_tree_' + gene + '_gene'
    output_file_name = 'consensus_tree_data_' + gene + '_gene'
    directory_name = gene + '_gene'
    intree_file_path = os.path.join(
        'output', directory_name, unrooted_tree_data_file_name)
    outtree_file_path = os.path.join('output', directory_name, outtree_file_name)
    output_file_path = os.path.join('output', directory_name, output_file_name)
    subprocess.call("./exec/consense")
    subprocess.call(["mv", "intree", intree_file_path])
    subprocess.call(["mv", "outtree", output_file_path])
    subprocess.call(["mv", "outfile", outtree_file_path])
