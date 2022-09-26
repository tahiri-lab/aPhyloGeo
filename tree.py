import subprocess
import pandas as pd
import os
import pathlib

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
file_name = '5seq/geo.csv'

specimen = 'id'   #"Please enter the name of the colum containing the specimens names: "

names = ['id', 'ALLSKY_SFC_SW_DWN', 'T2M', 'PRECTOTCORR', 'QV2M',
       'WS10M']


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

#-----------------------------------------

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

    # ecriture des matrices normalisees dans les fichiers respectifs
    with open(outfile_name, "w") as f:
        f.write("   " + str(len(nom_var)) + "\n")
        for j in range(nbr_seq):
            f.write(str(nom_var[j]))
            # petite boucle pour imprimer le bon nbr d'espaces
            for espace in range(11-len(str(nom_var[j]))):
                f.write(" ")
            for k in range(nbr_seq):
                # la normalisation se fait selon la formule suivante: (X - Xmin)/(Xmax - Xmin)
                f.write("{:.6f}".format(
                    (temp_tab[j][k] - min_value)/(max_value - min_value)) + " ")
            f.write("\n")
    subprocess.call(["rm", "outfile"])  # clean up

#-----------------------------------------

def create_tree(file_name, names):
    prepareDirectory()
    for i in range(1, len(names)):
        getDissimilaritiesMatrix(file_name, names[0], names[i], "infile") # liste a la position 0 contient les noms des specimens
        os.system("./exec/neighbor < input/input.txt")
        subprocess.call(["mv", "outtree", "intree"])
        subprocess.call(["rm", "infile", "outfile"])
        os.system("./exec/consense < input/input.txt" )
        newick_file = names[i].replace(" ", "_") + "_newick"
        subprocess.call(["rm", "outfile"])
        subprocess.call(["mv", "outtree", newick_file])
    subprocess.call(["rm", "intree"])




create_tree(file_name, names)

#prepareDirectory()
