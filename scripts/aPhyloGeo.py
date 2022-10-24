import subprocess   #! Deprecated import with BioPython integration
import pandas as pd
import os
import params #TODO : Migrate to .yaml format
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
import re
from csv import writer
import shutil
import Bio as Bio 
from Bio import SeqIO
from Bio import pairwise2

# ATTENTION AUX NOMS DES FICHIERS AVEC LES _
'''
bootstrap_threshold = 0
rf_threshold = 100
window_size = 5000
step_size = 500
data_names = ["T_min_à_2m_C_newick",
              "T_max_à_2m_C_newick"]
reference_gene_file = 'datasets/reference_gene.fasta'
'''
'''
bootstrap_threshold = 0
rf_threshold = 100
window_size = 10
step_size = 50
data_names = ['new_cases_smoothed_per_million_newick',
       'new_deaths_smoothed_per_million_newick', 'stringency_index_newick',
       'reproduction_rate_newick', 'people_vaccinated_per_hundred_newick',
       'people_fully_vaccinated_per_hundred_newick', 'population_density_newick',
       'median_age_newick', 'aged_65_older_newick', 'gdp_per_capita_newick',
       'cardiovasc_death_rate_newick', 'diabetes_prevalence_newick', 'female_smokers_newick',
       'male_smokers_newick', 'hospital_beds_per_thousand_newick'
    ]37
reference_gene_file = 'datasets/The_owid_final.fasta'
'''

bootstrap_threshold = params.bootstrap_threshold
rf_threshold = params.rf_threshold 
window_size = params.window_size
step_size = params.step_size
data_names = ['ALLSKY_SFC_SW_DWN_newick', 'T2M_newick', 'QV2M_newick', 'PRECTOTCORR_newick', 'WS10M_newick']

#reference_gene_file = 'datasets/The_37seq.fasta'

reference_gene_file = params.reference_gene_file


'''
file_name = 'donnees.csv'
specimen = 'Nom du specimen'   #"Please enter the name of the colum containing the specimens names: "
names = ['Nom du specimen','T min à 2m C',
        'T max à 2m C',
        'Humidité relative à 2m %']
'''


'''
file_name = 'Final_owid.csv'
specimen = 'Accession'   #"Please enter the name of the colum containing the specimens names: "
names = ['Accession','new_cases_smoothed_per_million',
       'new_deaths_smoothed_per_million', 'stringency_index',
       'reproduction_rate', 'people_vaccinated_per_hundred',
       'people_fully_vaccinated_per_hundred', 'population_density',
       'median_age', 'aged_65_older', 'gdp_per_capita',
       'cardiovasc_death_rate', 'diabetes_prevalence', 'female_smokers',
       'male_smokers', 'hospital_beds_per_thousand']
'''

file_name = params.file_name

specimen = params.specimen   #"Please enter the name of the colum containing the specimens names: "

names = params.names



def leastSquare(tree1, tree2):
    ls = 0.00
    # get all the terminal clades of the first tree (as example)
    specie_names = [specie.name for specie in tree1.get_terminals()]

    for i in range(tree1.count_terminals()-1):
        for j in range(i+1,tree1.count_terminals()):
            d1=(tree1.distance(tree1.find_any(specie_names[i]), tree1.find_any(specie_names[j])))
            d2=(tree2.distance(tree2.find_any(specie_names[i]), tree2.find_any(specie_names[j])))
            ls+=(abs(d1-d2))
    print(ls)
    return ls

def openCSV(nom_fichier_csv):
    df = pd.read_csv(nom_fichier_csv)
    return df

def getDissimilaritiesMatrix(df, column_with_specimen_name, column_to_search):
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

    # calculate des matrices normalisees 
    tab_df = pd.DataFrame(temp_tab)
    dm_df = (tab_df - min_value)/(max_value - min_value)
    dm_df = dm_df.round(6)

    matrix = [dm_df.iloc[i,:i+1].tolist() for i in range(len(dm_df))]
    dm = _DistanceMatrix(nom_var, matrix)
    return dm

def createTree(dm):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def climaticPipeline(file_name, names):
    trees = {}
    df = openCSV(file_name)
    for i in range(1, len(names)):
        dm = getDissimilaritiesMatrix(df, names[0], names[i]) # liste a la position 0 contient les noms des specimens
        trees[names[i]] = createTree(dm)

    leastSquare(trees[names[1]],trees[names[2]])

climaticPipeline(file_name, names)




def openFastaFile(reference_gene_file):
    sequences = {}
    with open(reference_gene_file) as sequencesFile:
        for sequence in SeqIO.parse(sequencesFile,"fasta"):
            sequences[sequence.id] = sequence.seq
    return sequences






def alignSequences(sequences):
        #print(sequences['ON129429'])
        #print(sequences['ON134852'])
    
    alignments = pairwise2.align.globalxx(sequences['ON129429'], sequences['ON134852'])
    print(alignments)
    #print(type(alignment))
    
    #subprocess.call(["./exec/muscle", "-in", reference_gene_file, "-physout", "infile", "-maxiters", "1", "-diags"])
    #f = open("infile", "r").read()
    #number_seq = int(f.split()[0])
    #subprocess.call(["cp", "infile", "alignment_result"])
    return 3



def createPhylogeneticTree(reference_gene_file, window_size, step_size, bootstrap_threshold, rf_threshold, data_names):
    #prepareDirectory()
    sequences = openFastaFile(reference_gene_file)
    print(sequences)
    number_seq = alignSequences(sequences)
    #slidingWindow(window_size, step_size)
    #files = os.listdir("output/windows")
    #for file in files:
    #    os.system("cp output/windows/" + file + " infile")
    #    createBoostrap()
    #    createDistanceMatrix()
    #    createUnrootedTree()
    #    createConsensusTree() # a modifier dans la fonction
    #    filterResults(reference_gene_file, bootstrap_threshold, rf_threshold, data_names, number_seq, file)
    

createPhylogeneticTree(reference_gene_file, window_size, step_size, bootstrap_threshold, rf_threshold, data_names)















#-----------------------------------------------------
#! Deprecated function with Biopython integration, to remove
def prepareDirectory():
    # delete the results of last analysis, if we have    ???
    with open("intree", "w"):
        pass
    # remove old newick files
    delete_path = os.listdir()

    for item in delete_path:
        if item.endswith("_newick"):
            os.remove(item)

    if os.path.exists("output/upload_gene.fasta") :
        os.remove("output/upload_gene.fasta")

#-----------------------------------------leaf







#-----------------------------------------  
def prepareDirectory():
    path_output_windows = './output/windows'                            
    isExist = os.path.exists(path_output_windows)

    if isExist:
        for f in os.listdir(path_output_windows):
            os.remove(os.path.join(path_output_windows, f))
    else:
        os.makedirs(path_output_windows)

    # delete the results of last analysis, if we have    
    delete_path = os.listdir('output')

    for item in delete_path:
        if item.endswith("_gene"):
            shutil.rmtree('output'+'/'+item)
        

    delete_path2 = os.listdir()

    for item in delete_path2:
        if item == "output.csv" or item.startswith("RAxML_") or item.startswith("outtree"):
            os.remove(item)
    
    with open('output.csv', 'w') as f:
        f.write("Gene,Arbre phylogeographique,Position ASM,Bootstrap moyen,RF normalise\n")

#prepareDirectory()

#--------------------------------------------------------------
#'1. Use the whole DNA sequences'

# main function


# number_seq = alignSequences(reference_gene_file)

# print(number_seq)

#-------------------------------------------------

def slidingWindow(window_size=0, step=0):
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
            # on connait la longueur de chaque sequence, 
            # donc on va recuperer chaque sequence et le retranscrire sur un autre fichier separes par un \n entre chaque
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
        # et on enregistre dans une list le nom en ordre des sequences
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
            with open("out", "r") as f, open("output/windows/" + str(debut) + "_" + str(fin), "w") as out:
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

def createBoostrap():
    filesize = os.path.getsize("infile")
    if filesize == 0:
        raise Exception("Infile for bootstrap was empty.")
    os.system("./exec/seqboot < input/bootstrap_input.txt")
    subprocess.call(["mv", "outfile", "infile"])

def createDistanceMatrix():
    os.system("./exec/dnadist < input/dnadist_input.txt")
    subprocess.call(["mv", "outfile", "infile"])

def createUnrootedTree():
    os.system("./exec/neighbor < input/neighbor_input.txt")
    subprocess.call(["rm", "infile", "outfile"])
    subprocess.call(["mv", "outtree", "intree"])


def createConsensusTree():
    os.system("./exec/consense < input/input.txt")
    # subprocess.call(["mv", "outtree", file])
    subprocess.call(["rm", "intree", "outfile"])

def filterResults(gene, bootstrap_threshold, rf_threshold, data_names, number_seq, aligned_file):
    bootstrap_average = calculateAverageBootstrap()
    if bootstrap_average < float(bootstrap_threshold):
        subprocess.call(["rm", "outtree"])
    else:
        for tree in data_names:
            #print(tree)
            calculateRfDistance(tree)
            rfn = standardizedRfDistance(number_seq)
            if rfn == None:                 #  '<=' not supported between instances of 'NoneType' and 'int'
                raise Exception(f'La distance RF n\'est pas calculable pour {aligned_file}.')                     # fix it 
            if rfn <= rf_threshold:
                runRaxML(aligned_file, gene, tree)
                cleanUp(aligned_file, tree)
                bootstrap_rax = calculateAverageBootstrapRax()
                if bootstrap_rax < float(bootstrap_threshold):
                    continue
                else:
                    calculateRfDistance(tree)
                    rfn_rax = standardizedRfDistance(number_seq)
                    if rfn_rax == None:     #  '<=' not supported between instances of 'NoneType' and 'int'
                        raise Exception(f'La distance RF pour Rax n\'est pas calculable pour {aligned_file}.')         # fix it 
                    if rfn_rax <= rf_threshold:
                        addToCsv(gene, tree, aligned_file, bootstrap_rax, rfn_rax)
                        keepFiles(gene, aligned_file, tree)
                        # a verifier ici
        subprocess.call(["rm", "outtree"])

def calculateAverageBootstrap():
    total = 0
    f = open("outtree", "r").read()
    numbers = re.findall(r'[)][:]\d+[.]\d+', f)
    for number in numbers:
        total = total + float(number[2:])
    if len(numbers) >0:
        average = total / len(numbers)
    else:
        average = -1
    return average

def calculateRfDistance(tree):
    os.system("cat " + tree + " >> infile")
    os.system("cat outtree >> infile")
    os.system("./exec/rf infile outfile tmp matrix")

def standardizedRfDistance(number_seq):
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

def runRaxML(aligned_file, gene, tree):
    current_dir = os.getcwd()
    file_name = os.path.basename(aligned_file + "_" + tree)
    input_path = os.path.join(current_dir, "output", "windows", aligned_file)

    # output_path = os.path.join(current_dir, "output", gene + "_gene")
    # IL FAUT CHANGER LE MODELE SELON LE GENE CHOISI
    os.system("./exec/raxmlHPC -s " + input_path + " -n " + file_name + " -N 100 -m GTRGAMMA -x 123 -f a -p 123")
    # output_path = os.path.join(output_path, file_name)
    # subprocess.call(["cp", input_path, output_path])

def cleanUp(file, tree):
    file = "RAxML_bipartitionsBranchLabels."+file+"_"+tree
    # directory = os.path.join("output", gene + "_gene", file)
    subprocess.call(["mv", file, "outtree"])
    files_to_delete = ['*bipartitions.*', '*bootstrap*', '*info*', '*bestTree*']
    for file in files_to_delete:
        os.system("rm -rf " +file)

def calculateAverageBootstrapRax():
    total = 0
    f = open("outtree", "r").read()
    numbers = re.findall(r'[\[]\d+[\]]', f)
    for number in numbers:
        total = total + float(number[1:(len(number)-1)])
    if len(numbers) > 0 :
        average = total / len(numbers)
    else:
        average = -1
    return average

def addToCsv(gene, tree, file, bootstrap_average, rfn):
    list = [gene, tree, file, bootstrap_average, rfn]
    with open('output.csv', 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow(list)
        f_object.close()

def keepFiles(gene, aligned_file, tree):
    current_dir = os.getcwd()
    file_name = os.path.basename(aligned_file + "_" + tree + "_tree")
    input_path = os.path.join(current_dir, "output", "windows", aligned_file)
    output_path = os.path.join(current_dir, "output", gene + "_gene")
    tree_path = os.path.join(output_path, file_name)
    subprocess.call(["cp", input_path, output_path]) # on garde l'ASM initial
    subprocess.call(["cp", "outtree", tree_path]) # on transfere l'arbre a garder dans le bon fichier
    subprocess.call(["mv", "output/windows/"+aligned_file+".reduced", output_path])




#displayGenesOption(window_size, step_size, bootstrap_threshold, rf_threshold, data_names,genes_chosen)