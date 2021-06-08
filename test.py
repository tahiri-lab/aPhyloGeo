from Bio import SeqIO
import subprocess
import re


def fetchingSequences():
    subprocess.run("./fetch_data.sh")

def changeNameSequences():
    sequences_file = open("output/sequences.fasta", "r")
    list_of_lines = sequences_file.readlines()
    for index in range(len(list_of_lines)):
        if list_of_lines[index].startswith(">"):
            splitted_line = list_of_lines[index].split("/")
            name = ">" + splitted_line[2] + "\n"
            list_of_lines[index] = name

    sequences_file = open("output/sequences.fasta", "w")
    sequences_file.writelines(list_of_lines)
    sequences_file.close()

def getORF1abGene():
    sequences_file = open("output/sequences.fasta", "r").read()
    list_of_sequences = sequences_file.split(">")
    s = 'ATCTAGGTTT(.*)TAACAACTAA'
    new_file = open("ORF1abGene.fasta", "w")
    for index in range(len(list_of_sequences)):
        if list_of_sequences[index] == "":
            continue
        name = list_of_sequences[index].split("\n")[0]
        gene_sequence = list_of_sequences[index].replace("\n", "")
        gene_sequence = (re.search(s, gene_sequence).group())
        new_file.writelines(">" + name + "\n")
        new_file.writelines(gene_sequence + "\n")
    
    new_file.close()


if __name__ == '__main__':
   fetchingSequences()
   changeNameSequences()
   getORF1abGene()
