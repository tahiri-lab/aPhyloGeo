import subprocess

def fetchingSequences():
    subprocess.run("./fetch_data.sh")

def changeNameSequences():
    sequences_file = open("sequences.fasta", "r")
    list_of_lines = sequences_file.readlines()
    for index in range(len(list_of_lines)):
        if list_of_lines[index].startswith(">"):
            splitted_line = list_of_lines[index].split("/")
            name = ">" + splitted_line[2] + "\n"
            list_of_lines[index] = name
    
    sequences_file = open("sequences.fasta", "w")
    sequences_file.writelines(list_of_lines)
    sequences_file.close()


def alignSequences():
    subprocess.run(["./muscle", "-in", "sequences.fasta", "-phyiout", "infile", "-maxiters", "1", "-diags"])

def createBoostrap():
    subprocess.run("./seqboot")
    subprocess.run(["mv", "outfile", "infile"])

def createDistanceMatrix():
    subprocess.run("./dnadist")


if __name__ == '__main__':
    # fetchingSequences()
    # changeNameSequences()
    # alignSequences()
    createBoostrap()
