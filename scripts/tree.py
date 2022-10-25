import pandas as pd
import os
import yaml
from yaml.loader import SafeLoader
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
import re
import toytree
import random
import toyplot.pdf

"""
Class description
"""

'''
file_name = 'donnees.csv'
specimen = 'Nom du specimen'   #"Please enter the name of the colum containing 
            the specimens names: "
names = ['Nom du specimen','T min à 2m C',
        'T max à 2m C',
        'Humidité relative à 2m %']
'''


'''
file_name = 'Final_owid.csv'
specimen = 'Accession'   #"Please enter the name of the colum containing the 
            specimens names: "
names = ['Accession','new_cases_smoothed_per_million',
       'new_deaths_smoothed_per_million', 'stringency_index',
       'reproduction_rate', 'people_vaccinated_per_hundred',
       'people_fully_vaccinated_per_hundred', 'population_density',
       'median_age', 'aged_65_older', 'gdp_per_capita',
       'cardiovasc_death_rate', 'diabetes_prevalence', 'female_smokers',
       'male_smokers', 'hospital_beds_per_thousand']
'''
#file_name = 'The_37_climate.csv'

# We open the params.yaml file and put it in the params variable
with open('params.yaml') as f:
    params = yaml.load(f, Loader=SafeLoader)
    print(params)
    
# Create variables from the yaml file content
file_name = params["file_name"]

specimen = params["specimen"]

names = params["names"]

#-------------------------------------------------------------------------------
#! Deprecated function, to remove
def prepareDirectory():
    # Delete the results of last analysis
    with open("intree", "w"):
        pass

    delete_path = os.listdir()

    for item in delete_path:
        if item.endswith("_newick"):
            os.remove(item)

    if os.path.exists("output/upload_gene.fasta") :
        os.remove("output/upload_gene.fasta")

#-------------------------------------------------------------------------------

def getDissimilaritiesMatrix(nom_fichier_csv, column_with_specimen_name, 
                            column_to_search):
    df = pd.read_csv(nom_fichier_csv)
    # Creation of a list containing the names of specimens and minimums 
    # tempratures
    meteo_data = df[column_to_search].tolist()
    nom_var = df[column_with_specimen_name].tolist()
    nbr_seq = len(nom_var)
    max_value = 0
    min_value = 0

    # First loop that allow us to calculate a matrix for each sequence
    temp_tab = []
    for e in range(nbr_seq):
        # A list that will contain every distances before normalisation
        temp_list = []
        for i in range(nbr_seq):
            maximum = max(float(meteo_data[e]), float(meteo_data[i]))
            minimum = min(float(meteo_data[e]), float(meteo_data[i]))
            distance = maximum - minimum
            temp_list.append(float("{:.6f}".format(distance)))

        # Allow to find the maximum and minimum value for the weather value and 
        # then to add the temporary list in an array   
        if max_value < max(temp_list):
            max_value = max(temp_list)
        if min_value > min(temp_list):
            min_value = min(temp_list)
        temp_tab.append(temp_list)

    # Calculate normalised matrix
    tab_df = pd.DataFrame(temp_tab)
    dm_df = (tab_df - min_value)/(max_value - min_value)
    dm_df = dm_df.round(6)

    matrix = [dm_df.iloc[i,:i+1].tolist() for i in range(len(dm_df))]
    dm = _DistanceMatrix(nom_var, matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    return tree

#-------------------------------------------------------------------------------

def leastSquare(tree1, tree2):
    """
    Method that calculates the least square distance between two trees.
    Trees must have the same number of leaves.
    Leaves must all have a twin in each tree.
    A tree must not have duplicate leaves
     x   x
    ╓╫╖ ╓╫╖
    123 312
 
    Args:
        tree1 (distanceTree object from biopython)
        tree2 (distanceTree object from biopython)
    
    Return:
        return result (double) the final distance between the two trees
    """
    result = 0.00
    leaves1 = tree1.get_terminals() #Produces a list of leaves from a tree
  
    leavesName = list(map(lambda l: l.name,leaves1))
 
    for i in leavesName:
        leavesName.pop(0)
        for j in leavesName:
            result1=(tree1.distance(tree1.find_any(i), tree1.find_any(j)))
            result2=(tree2.distance(tree2.find_any(i), tree2.find_any(j)))
            result+=(abs(result1-result2))
    print (result)
    return result

#-------------------------------------------------------------------------------

def draw_trees(trees):
    trees_newick= {}
    toytrees = []
    # Creating a multitree object from list of climatic trees
    for k,v in trees.items():
        trees_newick[k] = v.format('newick')
        ttree = toytree.tree(trees_newick[k], tree_format=1)
        toytrees.append(ttree)
    mtree = toytree.mtree(toytrees)

    # Setting up the stylings for nodes
    for tree in mtree.treelist:
        tree.style.edge_align_style={'stroke':'black','stroke-width':1}
        for node in tree.treenode.traverse():
            if node.is_leaf():
                node.add_feature('color', toytree.colors[7])  # terminals = grey
            else:
                node.add_feature('color', toytree.colors[1])  # internals/common = orange
    colors = tree.get_node_values('color', show_root=1, show_tips=1) 

    # Draw the climatic trees
    canvas, axes, mark = mtree.draw(nrows = round(len(mtree)/5), ncols=len(mtree), height=400, width=1000,node_sizes=8, node_colors=colors, tip_labels_align=True);
    for i in range(len(mtree)):
        rand_color = "#%03x" % random.randint(0, 0xFFF)
        axes[i].text(0,mtree.ntips,names[i+1],style={'fill':rand_color,'font-size':'10px', 'font-weight':'bold'});
    toyplot.pdf.render(canvas,'../viz/climactic_trees.pdf')
#-------------------------------------------------------------------------------

def create_tree(file_name, names):
    trees = {}
    for i in range(1, len(names)):
        trees[names[i]] = getDissimilaritiesMatrix(file_name, names[0], names[i])
        leaves = trees[names[i]].get_terminals()
        if i == 1:
            tree1 = trees[names[i]]
        if i == 2:
            tree2 = trees[names[i]]   
    leastSquare(trees[names[1]],trees[names[2]])
    draw_trees(trees)

#-------------------------------------------------------------------------------

create_tree(file_name, names)
