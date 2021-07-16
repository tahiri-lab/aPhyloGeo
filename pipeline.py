import subprocess
import os
import re
import pandas as pd
from scripts import reference

def menuGetTrees():
    names = []
    while True:
        count = input("How many climatic data tree will be used?: ")
        if not count.isnumeric():
            print("This is not a number.")
        elif int(count) < 1:
            print("The number cannot be lower than 1.")
        else:
            # VALIDER LE FORMAT NEWICK??
            try:
                for i in range(int(count)):
                    name = input("Name of the tree file (" + str(i+1) + "): " )
                    open(name, "r")
                    names.append(name)
            except:
                    print("This file does not exist or is empty.")
                    raise Exception
            break
    return names


def getBootstrapThreshold():
    valide = False
    while not valide:
        threshold = input("Enter the bootstrap value threshold between 0 and 100%: ")
        valide = threshold.isnumeric() and int(threshold) < 100 and int(threshold) > 0
        if not valide:
            print("Error, it must be a number between 1 and 100.")
        else:
            return threshold


def getRfThreshold():
    valide = False
    while not valide:
        threshold = input("Enter the Robinson and Foulds distance threshold between 0 and 100%: ")
        valide = threshold.isnumeric() and int(threshold) < 100 and int(threshold) > 0
        if not valide:
            print("Error, it must be a number between 1 and 100.")
        else:
            return threshold


def getSlidingWindowSize():
    while True:
        size = input("Sliding window size: ")
        try:
            size = int(size)
        except Exception:
            print("The window size must be a number.")
            continue
        if(size <= 0):
            print("The window size must be between greater than 0.")
        else:
            break
    return size

def getStepSize():
    while True:
        count = input("Step count: ")
        try:
            count = int(count)
        except Exception:
            print("The step count must be a number.")
            continue
        if(count <= 0):
            print("The window size must be between greater than 0.")
        else:
            break


def sliding_window(window_size=0, step=0):
    # Permet d'avoir le nombre de lignes totales dans le fichier
    f = open("infile", "r")
    line_count = -1
    for line in f:
        if line != "\n":
            line_count += 1
    f.close()
    f = open("infile", "r").read()
    # premier nombre de la premiere ligne du fichier represente le nbr de sequences
    num_seq = int((f.split("\n")[0]).split(" ")[0])
    # second nombre de la premiere ligne du fichier represente la longueur des sequences
    longueur = int((f.split("\n")[0]).split(" ")[1])
    # permet d'obtenir le nbr de lignes qui compose chaque sequence
    no_line = int(line_count/num_seq)

    # Recupere la sequence pour chaque variante
    with open("outfile", "w") as out:
        depart = 1
        fin = depart + no_line
        # on connait la longueur de chaque sequence, donc on va recuperer chaque sequence et le retranscrire sur un autre fichier separes par un \n entre chaque
        for i in range(0, int(num_seq)):
            f = open("infile", "r")
            lines_to_read = range(depart, fin)
            for position, line in enumerate(f):
                if position in lines_to_read:
                    out.write(line)
            out.write("\n")
            depart = fin
            fin = depart + no_line
    out.close()

    # on cree un fichier out qui contient chaque sequence sans espaces et on enregistre dans une list le nom en ordre des sequences
    with open("outfile", "r") as out, open("out", "w") as f:
        sequences = out.read().split("\n\n")
        list_names = []
        for seq in sequences:
            s = seq.replace("\n", " ").split(" ")
            if s[0] != "":
                list_names.append(s[0])
            s_line = s[1:len(seq)]
            for line in s_line:
                if line != "":
                    f.write(line)
            f.write("\n")
    out.close()
    f.close()

    # slide the window along the sequence
    debut = 0
    fin = debut + window_size
    while fin <= longueur:
        index = 0
        with open("out", "r") as f, open("output/windows/" + str(debut), "w") as out:
            out.write(str(num_seq) + " " + str(window_size) + "\n")
            for line in f:
                if line != "\n":
                    espece = list_names[index]
                    nbr_espaces = 11 - len(espece)
                    out.write(espece)
                    for i in range(nbr_espaces):
                        out.write(" ")
                    out.write(line[debut:fin] + "\n")
                    index = index + 1
        out.close()
        f.close()
        debut = debut + step
        fin = fin + step


def printOptionMenu():
    print('===============================================')
    print('Please select an option among the following: ')
    print('===============================================')
    print('1. Use the whole DNA sequences')
    print('2. Study specific genes of SARS-CoV-2')

def validateOptionMenu(window_size, step_size):
    while option != '1' or option != '2':
        option = input("Please enter 1 or 2: ")
        if option == '1':
            reference.getReferenceTree()
            break
        elif option == '2':
            print("ok2")
            break
        else:
            print('This is not a valid option.')


def menu(option=0):
    try:
        names = menuGetTrees()
        bootstrap_threshold = getBootstrapThreshold()
        rf_threshold = getRfThreshold()
        window_size = getSlidingWindowSize()
        step_size = getStepSize()
        printOptionMenu()
    except:
        return


def useWholeSequence(option=0):
    option = input("Use a sliding window?[y/n]\n")
    while True:
        if option == 'Y' or option == 'y':  # ici l'utilisateur choisit la fenetre coulissante
            size = input('Window size:\n')
            if int(size) > 0:
                print('Yes')
                break
        elif option == 'n' or option == 'N':  # ici, pas de fenetre coulissante on peut analyser toutes les sequences
            reference.get_reference_tree()
            break


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
    # subprocess.call(["mv", "infile", input_file_path])
    subprocess.call(["mv", "outfile", "infile"])


def getDissimilaritiesMatrix(nom_fichier_csv,column_with_specimen_name, column_to_search, outfile_name):
    df = pd.read_csv(nom_fichier_csv)
    # creation d'une liste contenant les noms des specimens et les temperatures min
    meteo_data = df[column_to_search].tolist()
    nom_var = df[column_with_specimen_name].tolist()
    nbr_seq = len(nom_var)
    # ces deux valeurs seront utiles pour la normalisation
    max_value = 0  
    min_value = 0

    # premiere boucle qui permet de calculer une matrice pour chaque sequence
    temp_tab = []
    for e in range(nbr_seq):
        # une liste qui va contenir toutes les distances avant normalisation
        temp_list = []
        for i in range(nbr_seq):        
            maximum = max(float(meteo_data[e]), float(meteo_data[i]))
            minimum = min(float(meteo_data[e]), float(meteo_data[i]))
            distance = maximum - minimum
            temp_list.append(float("{:.6f}".format(distance)))

        # permet de trouver la valeur maximale et minimale pour la donnee meteo et ensuite d'ajouter la liste temporaire a un tableau
        if max_value < max(temp_list):
            max_value = max(temp_list)
        if min_value > min(temp_list):
            min_value = min(temp_list)
        temp_tab.append(temp_list)
    
    # ecriture des matrices normalisees dans les fichiers respectifs
    with open(outfile_name, "w") as f:
        f.write("   " + str(len(nom_var)) + "\n")
        for j in range(nbr_seq):
            f.write(nom_var[j])
            # petite boucle pour imprimer le bon nbr d'espaces
            for espace in range(11-len(nom_var[j])):
                f.write(" ")
            for k in range(nbr_seq):
                # la normalisation se fait selon la formule suivante: (X - Xmin)/(Xmax - Xmin)
                f.write("{:.6f}".format((temp_tab[j][k] - min_value)/(max_value - min_value)) + " ")
            f.write("\n")
    subprocess.call(["rm", "outfile"]) # clean up


def createDistanceMatrix(gene):
    bootstrap_file_name = 'bootstrap_' + gene + '_gene'
    directory_name = gene + '_gene'
    input_file_path = os.path.join('output', directory_name, bootstrap_file_name)
    subprocess.call("./exec/dnadist")
    # subprocess.call(["cp", "infile", input_file_path])
    subprocess.call(["mv", "outfile", "infile"])

def createUnrootedTree(gene):
    distance_matrix_file_name = 'distance_matrix_' + gene + '_gene'
    directory_name = gene + '_gene'
    output_file_name = 'unrooted_tree_'+ gene + '_gene'
    input_file_path = os.path.join('output', directory_name, distance_matrix_file_name)
    output_file_path = os.path.join(
        'output', directory_name, output_file_name)
    subprocess.call("./exec/neighbor")
    # subprocess.call(["mv", "infile", input_file_path])
    subprocess.call(["mv", "outtree", "intree"])
    # subprocess.call(["mv", "outfile", output_file_path])


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
    # subprocess.call(["mv", "intree", intree_file_path])
    subprocess.call(["mv", "outtree", output_file_path])
    subprocess.call(["mv", "outfile", outtree_file_path])


if __name__ == '__main__':
    menu()
