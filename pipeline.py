from posix import listdir
import subprocess
import os
import re
import psutil
import pandas as pd
from csv import writer

def main():
    try:
        names = menu_get_trees()
        bootstrap_threshold = getBootstrapThreshold()
        rf_threshold = get_rf_threshold()
        window_size = get_sliding_window_size()
        step_size = get_step_size()
        validation_option_menu(window_size, step_size,
                           bootstrap_threshold, rf_threshold, names)
    except Exception as error:
        print(f'This error was caught: {error}.')
        # clean up all the windows even when an exception occured.
        # subprocess.call(["rm", "output/windows/*"])


def menu_get_trees():
    while True:
        count = input("How many climatic data tree will be used?: ")
        if not count.isnumeric():
            print("This is not a number.")
        elif int(count) < 1:
            print("The number cannot be lower than 1.")
        else:
            names = []
            for i in range(int(count)):
                name = input(f"Name of the tree file ({str(i+1)}): " )
                while name not in os.listdir("."):
                    print("This file does not exist or is empty.")
                    name = input(f"Name of the tree file ({str(i+1)}): " )
                names.append(name)
        return names


def getBootstrapThreshold():
    valide = False
    while not valide:
        threshold = input("Enter the bootstrap value threshold between 0 and 100%: ")
        valide = threshold.isnumeric() and int(threshold) <= 100 and int(threshold) >= 0
        if not valide:
            print("Error, it must be a number between 1 and 100.")
        else:
            return threshold


def get_rf_threshold():
    valide = False
    while not valide:
        threshold = input("Enter the Robinson and Foulds distance threshold between 0 and 100%: ")
        valide = threshold.isnumeric() and int(threshold) <= 100 and int(threshold) >= 0
        if not valide:
            print("Error, it must be a number between 1 and 100.")
        else:
            return int(threshold)


def get_sliding_window_size():
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

def get_step_size():
    while True:
        count = input("Step size: ")
        try:
            count = int(count)
        except Exception:
            print("The step size must be a number.")
            continue
        if(count <= 0):
            print("The step size must be between greater than 0.")
        else:
            return count


def sliding_window(window_size=0, step=0):
    # Permet d'avoir le nombre de lignes totales dans le fichier
    try:
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
            # on connait la longueur de chaque sequence, donc on va recuperer
            #  chaque sequence et le retranscrire sur un autre fichier separes
            #  par un \n entre chaque
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

        # on cree un fichier out qui contient chaque sequence sans espaces
        #  et on enregistre dans une list le nom en ordre des sequences
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
            with open("out", "r") as f, open(f"output/windows/{str(debut)}_{str(fin)}", "w") as out:
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
    except:
        print("An error occurred.")

    # clean up
    os.system("rm out outfile infile")


def validation_option_menu(window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    print('''
    ===============================================
    Please select an option among the following:
    ===============================================
    1. Use the whole DNA sequences
    2. Study specific genes of SARS-CoV-2'''
    )
    while True:
        option = input("Please enter 1 or 2: ")
        if option == '1':
            gene = 'reference'
            create_phylo_tree(gene, window_size, step_size, bootstrap_threshold, rf_threshold, data_names)
            break
        elif option == '2':
            display_gene_options(window_size, step_size,
                               bootstrap_threshold, rf_threshold, data_names)
            break
        else:
            print('This is not a valid option.')


def display_gene_options(window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    options = {
             "1": "ORF1ab","2": "S",
             "3": "ORF3a", "4": "ORF3b",
             "5": "E", "6": "M",
             "7": "ORF6", "8": "ORF7a",
             "9": "ORF7b", "10": "ORF8",
             "11": "N", "12": "ORF10"
             }

    print(
        '''================================================================================
           Choose among the following genes to analyze seperated by spaces (ex: 1 8 11): 
           ================================================================================''')
    for number, gene in options.items():
        print(number, ":", gene)
    while True:
        values = input()
        values = values.split()
        # d'abord verifier que toutes les valeurs sont valides
        valides = False
        for value in values:
            if options.get(value) == None:
                print("This is not a valid input. Example of usage: 1 5 12. ")
                valides = False
                break
            else:
                valides = True
        if valides:
            genes = {
                    'ORF1ab': 'ATGGAGAGCC(.*)TAACAACTAA',
                    'S': 'ATGTTTGTTT(.*)TTACACATAA',
                    'ORF3a': 'ATGGATTTGT(.*)GCCTTTGTAA',
                    'ORF3b': 'ATGAGGCTTT(.*)GCCTTTGTAA',
                    'E': 'ATGTACTCAT(.*)TCTGGTCTAA',
                    'M': 'ATG[GT]CAGATT(.*)TGTACAGTAA',
                    'ORF6': 'ATGTTTCATC(.*)GATTGA[CT]TAA',
                    'ORF7a': 'ATGAAAATTAT(.*)GACAGAATGA',
                    'ORF7b': 'ATGATTGAACTTTCATTAATTGACTTCTATTTGTGCTTTTTAGCCTTTCTGCTATTCCTTGTTTTAATTATGCTTATTATCTTTTGGTTCTCACTTGAACTGCAAGATCATAATGAAACTTGTCACGCCTAA',
                    'ORF8': 'ATGAAATTTCTTGTTTT(.*)TTT[TC]ATCTAA',
                    'N': 'ATGTCT[CG][AT][TA]AAT(.*)TCAGGCCTAA',
                    'ORF10': 'ATGGGCTATA(.*)TCTCACATAG'}
            for value in values:
                gene = options.get(value)
                pattern = genes.get(gene)
                get_gene(gene, pattern)
                create_phylo_tree(gene, window_size, step_size, bootstrap_threshold, rf_threshold, data_names)
            subprocess.call(["make", "clean"])
            break


def get_gene(gene, pattern): 
    sequences_file = open("output/reference_gene.fasta", "r").read()
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


def create_phylo_tree(gene, window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    number_seq = align_sequence(gene)
    sliding_window(window_size, step_size)
    files = os.listdir("output/windows")
    for file in files:
        try:
            os.system("cp output/windows/" + file + " infile")
            create_bootstrap()
            run_dnadist()
            run_neighbor()
            run_consense() # a modifier dans la fonction
            filter_results(gene, bootstrap_threshold, rf_threshold, data_names, number_seq, file)
            subprocess.call(["rm", "output/windows/"+file])
        except Exception as error:
            raise error

def align_sequence(gene):
    sequences_file_name = gene + '_gene.fasta'
    directory_name = gene + '_gene'
    if gene == 'reference':
        file_path = os.path.join('output', sequences_file_name)
    else:
        file_path = os.path.join('output', directory_name, sequences_file_name)
    subprocess.call(["./exec/muscle", "-in", file_path, "-physout", "infile", "-maxiters", "1", "-diags"])
    file_path =os.path.join('output', directory_name)
    # subprocess.call(["cp", "infile", file_path])
    f = open("infile", "r").read()
    number_seq = int(f.split()[0])
    return number_seq


def create_bootstrap():
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Infile for bootstrap was empty.")
    os.system("./exec/seqboot < input/bootstrap_input.txt")
    subprocess.call(["mv", "outfile", "infile"])

    # check if the bootstrap file is empty
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Bootstrap file is empty.")


def create_matrix(nom_fichier_csv,column_with_specimen_name, column_to_search, outfile_name):
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

        # permet de trouver la valeur maximale et minimale pour la donnee meteo
        #  et ensuite d'ajouter la liste temporaire a un tableau
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


def run_dnadist():
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Infile for distance was empty.")
    os.system("./exec/dnadist < input/dnadist_input.txt")
    subprocess.call(["mv", "outfile", "infile"])
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Something went wrong with distance file.")

    


def run_neighbor():
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Infile for unrooted tree was empty.")
    os.system("./exec/neighbor < input/neighbor_input.txt")
    subprocess.call(["rm", "infile", "outfile"])
    subprocess.call(["mv", "outtree", "intree"])
    filesize = os.path.getsize("intree")
    if filesize == 0:
        raise Exception("Something went wrong with neighbor file.")


def run_consense():
    filesize = os.path.getsize("intree")
    if filesize == 0:
        raise Exception("Intree for consense was empty.")
    os.system("./exec/consense < input/input.txt")
    # subprocess.call(["mv", "outtree", file])
    subprocess.call(["rm", "intree", "outfile"])


def calculate_average_bootstrap():
    total = 0
    f = open("outtree", "r").read()
    numbers = re.findall(r'[)][:]\d+[.]\d+', f)
    for number in numbers:
        total = total + float(number[2:])
    average = total / len(numbers)
    return average


def calculate_average_bootstrap_rax():
    total = 0
    f = open("outtree", "r").read()
    numbers = re.findall(r'[\[]\d+[\]]', f)
    for number in numbers:
        total = total + float(number[1:(len(number)-1)])
    average = total / len(numbers)
    return average


def calculate_rf_distance(tree):
    filesize = os.path.getsize('outtree')
    if filesize == 0:
        raise Exception('Outree pour le calcul de distance est vide.')
    os.system(f"cat {tree} >> infile")
    os.system("cat outtree >> infile")
    os.system("./exec/rf infile outfile tmp matrix")

def standardized_rf_distance(number_seq):
    # clean up the repository
    subprocess.call(["rm", "infile", "matrix", "tmp"])
    # find the rf
    f = open("outfile", "r").read()
    words = re.split(r'[ \n]', f)
    for i in range(len(words)):
        if words[i] == "=":
            rf = int(words[i+1])
            normalized_rf = (rf/(2*number_seq-6))*100
            subprocess.call(["rm", "outfile"])
            return normalized_rf


def run_raxml(aligned_file, tree):
    current_dir = os.getcwd()
    file_name = os.path.basename(f"{aligned_file}_{tree}")
    input_path = os.path.join(current_dir, "output", "windows", aligned_file)
    # trouver le nombre de cpu de la machine
    nbr_cpu = psutil.cpu_count(logical = False)
    # n'a pas pu trouver le nbr de cpu, on roule la version sequentielle de raxML
    if nbr_cpu == None: 
        os.system(
            f"./exec/raxmlHPC -s {input_path}  -n {file_name}  -N 100 -m GTRGAMMA -x 123 -f a -p 123")
    else: # on roule la version PTHREAD en specifiant le nombre de cpu
        os.system(
            f"./exec/raxmlHPC-PTHREADS -s {input_path} -n {file_name} -N 100 -m GTRGAMMA -x 123 -f a -p 123 -T {nbr_cpu}")
    # output_path = os.path.join(output_path, file_name)
    # subprocess.call(["cp", input_path, output_path])


def filter_results(gene, bootstrap_threshold, rf_threshold, data_names, number_seq, aligned_file):
        bootstrap_average = calculate_average_bootstrap()
        if bootstrap_average < float(bootstrap_threshold):
            subprocess.call(["rm", "outtree"])
        else:
            for tree in data_names:
                calculate_rf_distance(tree)
                rfn = standardized_rf_distance(number_seq)
                if rfn == None:      #'<=' not supported between instances of 'NoneType' and 'int'
                    raise Exception(f'La distance RF n\'est pas calculable pour {aligned_file}.')     #fix it 
                if rfn <= rf_threshold:
                    run_raxml(aligned_file, tree)
                    clean_up(aligned_file, tree)
                    bootstrap_rax = calculate_average_bootstrap_rax()
                    if bootstrap_rax < float(bootstrap_threshold):
                        continue
                    else:
                        calculate_rf_distance(tree)
                        rfn_rax = standardized_rf_distance(number_seq)
                        if rfn_rax == None:         #'<=' not supported between instances of 'NoneType' and 'int'
                            raise Exception(f'La distance RF pour Rax n\'est pas calculable pour {aligned_file}.')  # fix it
                        if rfn_rax <= rf_threshold:
                            add_to_csv(gene, tree, aligned_file, bootstrap_rax, rfn_rax)
                            keep_files(gene, aligned_file, tree)
                            # a verifier ici
            subprocess.call(["rm", "outtree"])


def keep_files(gene, aligned_file, tree):
    current_dir = os.getcwd()
    file_name = os.path.basename(f"{aligned_file}_{tree}_tree")
    input_path = os.path.join(current_dir, "output", "windows", aligned_file)
    output_path = os.path.join(current_dir, "output", gene + "_gene")
    tree_path = os.path.join(output_path, file_name)
    subprocess.call(["cp", input_path, output_path]) # on garde l'ASM initial
    subprocess.call(["cp", "outtree", tree_path]) # on transfere l'arbre a garder dans le bon fichier


def add_to_csv(gene, tree, file, bootstrap_average, rfn):
    list = [gene, tree, file, bootstrap_average, rfn]
    with open('output.csv', 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow(list)
        f_object.close()


def clean_up(file, tree):
    reduced_file = f"{file}.reduced"
    file = f"RAxML_bipartitionsBranchLabels.{file}_{tree}"
    # directory = os.path.join("output", gene + "_gene", file)
    subprocess.call(["mv", file, "outtree"])
    files_to_delete = ['*bipartitions.*',
                       '*bootstrap*',
                       '*bestTree*',
                       '*info*']
    for file in files_to_delete:
        os.system("rm -rf " +file)
    subprocess.call(["rm", f"output/windows/{reduced_file}"])

if __name__ == '__main__':
    main()
