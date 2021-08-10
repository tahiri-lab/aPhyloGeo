# Phylotree

Cet outil est un pipeline bioinformatique dédié à l'analyse du SARS-CoV-2. Dans un premier temps, il permet d'obtenir des arbres phylogéographiques à partir des données météorologiques des régions où les échantillons du SARS-CoV-2 ont été récoltés. Dans un deuxième temps, cet outil permet de comparer des arbres phylogénétiques obtenus à partir d'alignements de séquences multiples (ASM) des variantes du SARS-CoV-2 avec les arbres phylogéographiques afin de conserver les régions des séquences qui présentent une ressemblance topologique avec l'arbre phylogéographique, et donc, qui pourrait potentiellement être corrélé avec cette variable météorologique.

## Workflow
Pour qu

## Comment l'utiliser?

### Pré-requis
Avant d'utiliser ce programme, assurez-vous d'avoir installé toutes les librairies nécessaires à son bon fonctionnement. Pour se faire, simplement taper la commande suivante :

```
pip install -r requirements.txt
```
### Création des arbres phylogéographiques

Afin d'obtenir les arbres phylogéographiques, il faut créer un fichier csv contenant les données d'intérêt (voir l'exemple [donnees.csv](./data/donnees.csv))

Ce fichier devrait avoir une structure similaire au fichier exemple. Une colonne contenant les noms des spécimens devrait être présente puisque le nom de cette colonne sera demandé au lancement du programme.

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
Please enter the name of the csv file (this is a relative path): data
Could not find the csv file with this name or the file is empty.
Please enter the name of the csv file: data/donnees.csv
Number of trees to create: 
5
Please enter the name of the colum containing the specimens names: Nom du specimens
Error, this column name does not exist.
Please enter the name of the colum containing the specimens names: Nom du specimen
Please enter the name of the column to analyze in your csv file (1): T min à 2m (C)
Please enter the name of the column to analyze in your csv file (2): T max à 2m (C)
Please enter the name of the column to analyze in your csv file (3): Humidité relative à 2m(%)
Please enter the name of the column to analyze in your csv file (4): Précipitation totale sur le mois (mm)
Please enter the name of the column to analyze in your csv file (5): Pression en surface (kPa)
```

```Assurez-vous que les noms des colonnes soient identiques à celles que vous avez saisies pour que le programme puisse les retrouvés! Les noms des séquences à étudier ainsi que les noms des spécimens dans ce fichier doivent également correspondre exactement.```

### Création des arbres phylogénétiques et analyses phylogéographiques

<p>
1. Les fichiers des séquences à analyser (.fasta) doivent être mis dans le dossier [data](./data). À noter que l'entête des fichiers de séquences doivent respecter ce format :

```>hCoV-19/Nigeria/S38/2020|EPI_ISL_2399462|2020-10-12```

Afin de pouvoir différencier chaque séquence avec un nom différent, cet outil va isoler le 3e mot rencontré dans l'entête. Dans cet exemple, la séquence sera alors nommée ```S38```. Dans le cas où ce mot serait plus long que 10 caractères, le programme ne prendra que les 10 premiers. Par exemple, si nous avons ceci comme entête:

```>hCoV-19/India/GJ-GBRC560b/2021|EPI_ISL_1677798|2021-01-12```

Dans ce cas, le nom de la séquence sera seulement ```GJ-GBRC560``` puisqu'il est composé d'exactement 10 caractères. Ceci est important à comprendre, surtout pour les noms que vous donnerez aux séquences dans le fichier ```csv``` qui doivent être identiques à ces noms en question.

</p>




