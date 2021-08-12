# Phylotree

```Phylotree``` is a bioinformatics pipeline dedicated to the analysis of ```SARS-CoV-2```. On the one hand, it allows to obtain trees from meteorological data of the regions where ```SARS-CoV-2``` samples have been collected. On the other hand, this tool allows topological comparison of phylogenetic trees from multiple sequence alignments (MSAs) with trees from weather data using the [Robinson-Foulds (RF) metric](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub). ASMs that yield trees with a significant ```RF``` value are stored in folders with the alignment and tree in question.

## Workflow

![](./img/workflow_en.jpeg)

## How to use?

### Prerequisites
Before using this program, make sure that you have installed all the necessary libraries for it to work properly. To do this, simply type the following command:

```
pip install -r requirements.txt
```

### Creation of phylogeographic trees

In order to obtain the phylogeographic trees, you need to create a csv file containing the data of interest (see the example [data.csv](./data/data.csv))

This file should have a structure similar to the example file. A column containing the names of the specimens should be present since the name of this column will be requested when the program is launched. **The column names should not contain parentheses!

Then the other columns should contain the meteorological variables that are to be studied. The values must be **numeric** for the program to work properly.

Once the file is created with the data, it is now time to create the trees with this command: 

```
make tree
```
Here is an example of what is displayed on the terminal and the entries entered according to the example file [data.csv](./data/data.csv) :

```
====================================================================================================================
Before running this script, please make sure there is a .csv file containing all the data to analyze in this repo
====================================================================================================================
Please enter the name of the csv file (this is a relative path): data/donnees.csv
Number of trees to create: 
5
Please enter the name of the colum containing the specimens names: Nom du specimen
Please enter the name of the column to analyze in your csv file (1): T min à 2m C
Please enter the name of the column to analyze in your csv file (2): T max à 2m C
Please enter the name of the column to analyze in your csv file (3): Humidité relative à 2m %
Please enter the name of the column to analyze in your csv file (4): Précipitation totale sur le mois mm
Please enter the name of the column to analyze in your csv file (5): Pression en surface kpa
```

Make sure that the column names are identical to the ones you entered so that the program can find them! The names of the sequences to be studied and the names of the specimens in this file must also match exactly.

The phylogeographic trees should then be created in the current folder with the column names followed by ``newick`` (e.g. ``T_max_to_2m_C_newick``). 


### Creating phylogenetic trees and phylogeographic analyses

1. The files of the sequences to be analyzed (.fasta) must be put in the folder [data](./data). Note that the header of the sequence files must respect this format (the sequences obtained on [GISAID](https://www.gisaid.org) all respect this format):

```>hCoV-19/Nigeria/S38/2020|EPI_ISL_2399462|2020-10-12```

In order to be able to differentiate each sequence with a different name, this tool will isolate the 3rd word encountered in the header. In this example, the sequence will be named "S38". In case this word is longer than 10 characters, the program will only take the first 10. For example, if we have this as a header:

```>hCoV-19/India/GJ-GBRC560b/2021|EPI_ISL_1677798|2021-01-12```

In this case, the sequence name will only be ``GJ-GBRC560`` since it is exactly 10 characters long. This is important to understand, especially for the names you will give to the sequences in the ``csv`` file, which must be identical to the single names in question.

2. Once the sequences have been placed in the [data] folder (./data), all you have to do is run the program with the command ``make''. Here is an example of what should appear on the terminal and an example of the entries:

```
How many climatic data tree will be used?: 1
Name of the tree file (1): T_max_à_2m_C_newick
Enter the bootstrap value threshold between 0 and 100%: 10
Enter the Robinson and Foulds distance threshold between 0 and 100%: 100
Sliding window size: 100
Step count: 10
===============================================
Please select an option among the following: 
===============================================
1. Use the whole DNA sequences
2. Study specific genes of SARS-CoV-2
Please enter 1 or 2: 2
================================================================================
Choose among the following genes to analyze seperated by spaces (ex: 1 8 11): 
================================================================================
1 : ORF1ab
2 : S
3 : ORF3a
4 : ORF3b
5 : E
6 : M
7 : ORF6
8 : ORF7a
9 : ORF7b
10 : ORF8
11 : N
12 : ORF10
7
```
Once the program is run, depending on the bootstrap value and the threshold distance, the phylogenetic trees of interest with their corresponding ASM (multiple sequence alignment) will be kept in the folders of each selected [gene](./output) or in the [reference] folder (./output/reference_gene) if the complete sequence has been studied. Finally, the [output.csv](output.csv) file will contain all the significant ASMs selected with their gene in question, their position on the sequence, their bootstrap value and its rf distance.


## Potential problems encountered

+ For `macOS` users, it is likely that your computer is blocking access to the `MUSCLE` program. If this is the case, simply go to the privacy settings on your machine and give the program access

## References

1. Sequence alignment tool : `MUSCLE`
    + [Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797.](https://academic.oup.com/nar/article/32/5/1792/2380623)
    doi:10.1093/nar/gkh340]
    + [Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity BMC Bioinformatics, (5) 113.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113)
    doi:10.1186/1471-2105-5-113

2. Phylip Package: `Seqboot`, `DNADist`, `Neighbor`, `Consense`
    + [J. Felsenstein. 1989. PHYLIP - Phylogeny Inference Package (Version 3.2) . Cladistics. 5: 164-166.](https://evolution.genetics.washington.edu/phylip.html)

3. Robinson-Foulds metric:
    + [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)

4. Phylogenetic analysis: `RAxML`
    + [A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053?login=true)


------

```Phylotree``` est un pipeline bioinformatique dédié à l'analyse du ```SARS-CoV-2```. D'une part, il permet d'obtenir des arbres issus de données métérologiques des régions où les échantillons du SARS-CoV-2 ont été récoltés. D'autre part, cet outil permet de comparer topologiquement des arbres phylogénétiques issus des alignements de séquences multiples (ASM) aux arbres issus des données météo à l'aide du calcul de la [distance de Robinson et Foulds (RF)](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub). Les ASM qui donnent des arbres avec une valeur ```RF``` significative sont conservés dans des dossiers avec l'alignement et l'arbre en question.

## Workflow
![](./img/workflow_fr.jpeg)

## Comment l'utiliser?

### Pré-requis
Avant d'utiliser ce programme, assurez-vous d'avoir installé toutes les librairies nécessaires à son bon fonctionnement. Pour se faire, simplement taper la commande suivante :

```
pip install -r requirements.txt
```
### Création des arbres phylogéographiques

Afin d'obtenir les arbres phylogéographiques, il faut créer un fichier csv contenant les données d'intérêt (voir l'exemple [donnees.csv](./data/donnees.csv))

Ce fichier devrait avoir une structure similaire au fichier exemple. Une colonne contenant les noms des spécimens devrait être présente puisque le nom de cette colonne sera demandé au lancement du programme. **Les noms des colonnes ne doivent pas contenir de parenthèses!**

Ensuite, les autres colonnes devraient contenir les variables métérologiques qui sont à étudier. Les valeurs doivent être **numériques** pour que le programme puisse fonctionner correctement.

Une fois le fichier créé avec les données, il est maintenant temps de créer les arbres avec cette commande : 

```
make tree
```
Voici un exemple de ce qui s'affiche sur le terminal et les saisies entrées selon le fichier exemple [donnees.csv](./data/donnees.csv) :
```
====================================================================================================================
Before running this script, please make sure there is a .csv file containing all the data to analyze in this repo
====================================================================================================================
Please enter the name of the csv file (this is a relative path): data/donnees.csv
Number of trees to create: 
5
Please enter the name of the colum containing the specimens names: Nom du specimen
Please enter the name of the column to analyze in your csv file (1): T min à 2m C
Please enter the name of the column to analyze in your csv file (2): T max à 2m C
Please enter the name of the column to analyze in your csv file (3): Humidité relative à 2m %
Please enter the name of the column to analyze in your csv file (4): Précipitation totale sur le mois mm
Please enter the name of the column to analyze in your csv file (5): Pression en surface kpa
```

```Assurez-vous que les noms des colonnes soient identiques à celles que vous avez saisies pour que le programme puisse les retrouvés! Les noms des séquences à étudier ainsi que les noms des spécimens dans ce fichier doivent également correspondre exactement.```

Les arbres phylogéographiques devraient alors être créés dans le dossier courant avec le nom des colonnes suivi de ```_newick``` (ex: ```T_max_à_2m_C_newick```). 


### Création des arbres phylogénétiques et analyses phylogéographiques

1. Les fichiers des séquences à analyser (.fasta) doivent être mis dans le dossier [data](./data). À noter que l'entête des fichiers de séquences doivent respecter ce format ( les séquences obtenues sur [GISAID](https://www.gisaid.org) respectent tous ce format ) :

```>hCoV-19/Nigeria/S38/2020|EPI_ISL_2399462|2020-10-12```

Afin de pouvoir différencier chaque séquence avec un nom différent, cet outil va isoler le 3e mot rencontré dans l'entête. Dans cet exemple, la séquence sera alors nommée ```S38```. Dans le cas où ce mot serait plus long que 10 caractères, le programme ne prendra que les 10 premiers. Par exemple, si nous avons ceci comme entête:

```>hCoV-19/India/GJ-GBRC560b/2021|EPI_ISL_1677798|2021-01-12```

Dans ce cas, le nom de la séquence sera seulement ```GJ-GBRC560``` puisqu'il est composé d'exactement 10 caractères. Ceci est important à comprendre, surtout pour les noms que vous donnerez aux séquences dans le fichier ```csv``` qui doivent être identiques à ces noms isolés en question.

2. Une fois les séquences bien placées dans le dossier [data](./data), il suffit maintenant de lancer le programme avec la commande ```make```. Voici un exemple de ce qui devrait apparaître sur le terminal ainsi qu'un exemple des saisies : 

```
How many climatic data tree will be used?: 1
Name of the tree file (1): T_max_à_2m_C_newick
Enter the bootstrap value threshold between 0 and 100%: 10
Enter the Robinson and Foulds distance threshold between 0 and 100%: 100
Sliding window size: 100
Step count: 10
===============================================
Please select an option among the following: 
===============================================
1. Use the whole DNA sequences
2. Study specific genes of SARS-CoV-2
Please enter 1 or 2: 2
================================================================================
Choose among the following genes to analyze seperated by spaces (ex: 1 8 11): 
================================================================================
1 : ORF1ab
2 : S
3 : ORF3a
4 : ORF3b
5 : E
6 : M
7 : ORF6
8 : ORF7a
9 : ORF7b
10 : ORF8
11 : N
12 : ORF10
7
```
Une fois le programme lancé, selon la valeur de ```bootstrap``` et la distance ```rf``` seuil, les arbres phylogénétiques d'intérêts avec leur ASM (alignement de séquence multiple) correspondant seront gardés dans les dossiers de chaque [gène](./output) sélectionnés ou dans le dossier [référence](./output/reference_gene) si la séquence complète a été étudiée. Finalement, le fichier [output.csv](output.csv) contiendra tous les ASM significatifs retenus avec leur gène en question, leur position sur la séquence, leur valeur de ```bootstrap``` et sa distance ```rf```.


## Potentiels problèmes rencontrés

+ Pour les utilisateurs de ```macOS```, il est probable que votre ordinateur bloque l'accès au programme `MUSCLE`. Si c'est le cas, simplement aller dans les paramètres de confidentialité de votre machine et donner l'accès au programme.

## Références
Plusieurs programmes externes ont été utilisé dans l'élaboration de ce pipeline :

1. Programme pour l'alignement des séquences: ```MUSCLE``` 
    + [Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797.](https://academic.oup.com/nar/article/32/5/1792/2380623)
    doi:10.1093/nar/gkh340]
    + [Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity BMC Bioinformatics, (5) 113.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113)
    doi:10.1186/1471-2105-5-113

2. Paquet Phylip: ```Seqboot```, ```DNADist```, ```Neighbor```, ```Consense```
    + [J. Felsenstein. 1989. PHYLIP - Phylogeny Inference Package (Version 3.2) . Cladistics. 5: 164-166.](https://evolution.genetics.washington.edu/phylip.html)

3. Calcul de la distance de ```Robinson et Foulds``` : [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)

4. Outil pour les analyses phylogénétiques : ```RAxML```
    + [A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053?login=true)


