import subprocess
import pandas as pd
import os
import params
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
import re
#from ete3 import Tree
#Bonjour Frank!

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
#file_name = 'The_37_climate.csv'
file_name = params.file_name

specimen = params.specimen   #"Please enter the name of the colum containing the specimens names: "

names = params.names


#-----------------------------------------------------
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

def getDissimilaritiesMatrix(nom_fichier_csv, column_with_specimen_name, column_to_search, outfile_name):
    df = pd.read_csv('datasets/' + nom_fichier_csv)
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
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    ##print(tree)
    #Phylo.write(tree, outfile_name, "newick")
    #tree = Phylo.read(outfile_name, "newick")
    #print(type(tree))
    return tree

#-----------------------------------------

def leastSquare(tree1, tree2):
    result = 0.00
    leaves1 = tree1.get_terminals()
    
    leavesName = []
    for leave in leaves1:
        print(leaves1)
        leavesName.append(leave.name)

    print(leavesName)
    leavesNameTemp = leavesName

    for i in leavesName:
        leavesNameTemp.pop(0)
        print("------------------")
        for j in leavesNameTemp:
            result1=(tree1.distance(tree1.find_any(i), tree1.find_any(j)))
            result2=(tree2.distance(tree2.find_any(i), tree2.find_any(j)))
            print(abs(result1-result2))
            result+=(abs(result1-result2))
            
    print("************************")

    print(result)
    return result


def create_tree(file_name, names):
    #prepareDirectory()
    trees = {}
    for i in range(1, len(names)):
        trees[names[i]] = getDissimilaritiesMatrix(file_name, names[0], names[i], "infile") # liste a la position 0 contient les noms des specimens
        ##print(trees[names[i]].count_terminals())
        leaves = trees[names[i]].get_terminals()
        ##print(leaves[0])
        ##print(leaves[1])
        ##print(trees[names[i]].distance(leaves[0], leaves[1]))
        #print(trees[names[i]])
        if i == 1:
            tree1 = trees[names[i]]
        if i == 2:
            tree2 = trees[names[i]]   

        #os.system("./exec/neighbor < input/input.txt")
        #subprocess.call(["mv", "outtree", "intree"])
        #subprocess.call(["rm", "infile", "outfile"])
        #os.system("./exec/consense < input/input.txt" )
        #newick_file = names[i].replace(" ", "_") + "_newick"
        #subprocess.call(["rm", "outfile"])
        #subprocess.call(["mv", "outtree", newick_file])
    #subprocess.call(["rm", "intree"])

    leastSquare(trees[names[1]],trees[names[2]])
    #print(tree1.robinson_foulds(tree2))

create_tree(file_name, names)

#prepareDirectory()