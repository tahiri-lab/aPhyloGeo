﻿﻿﻿﻿﻿<h1  align="center"> aPhyloGeo <p align='center'>
        [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
        [![Contributions](https://img.shields.io/badge/contributions-welcome-blue.svg)](https://pysd.readthedocs.io/en/latest/development/development_index.html)
        [![Py version](https://img.shields.io/pypi/pyversions/pysd.svg)](https://pypi.python.org/pypi/pysd/)
        [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftahiri-lab%2FaPhylogeo&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
        [![GitHub release](https://img.shields.io/github/v/release/tahiri-lab/aPhylogeo.svg?maxAge=3600)](https://github.com/tahiri-lab/aPhylogeo/releases/)
        </p>


<h2  align="center"> 🌳 Application multiplateforme pour l'analyse des arbres phylogénétiques avec des paramètres climatiques</h2>

<details open>
  <summary>Table des matières</summary>
  <ol>
    <li>
      <a href="#about-the-project">À propos du projet</a>
    </li>
    <li>
      <a href="#Installation">Installation</a>
      <ul>
        <li><a href="#Linux-UNIX-and-Mac-OS-versions">Linux/UNIX and Mac OS versions</a></li>
      </ul>
    </li>
    <!--<li> Analyses disponibles</li>
      <ul>
        <li><a href="#Group-creation">Création de groupe</a></li>
        <li><a href="#SimPlot-analysis">Analyse SimPlot</a></li>
        <li><a href="#Similarity-networks">Réseaux de similarité</a></li>
        <li><a href="#BootScan-analysis">Analyse BootScan</a></li>
        <li><a href="#Findsites">Findsites</a></li>
        <li><a href="#Detection-of-recombination">Détection de recombinaison</a></li>
      </ul>-->
     <li>
      <a href="#Settings">Paramètres</a>
    </li>
    <li>
      <a href="#Example">Exemple</a>
      <ul>
        <li><a href="#Input">Entrée</a></li>
        <li><a href="#Output">Sortie</a></li>
      </ul>
    </li>
    <li>
      <a href="#References">Références</a>
    </li>
    <li>
      <a href="#contact">Contact</a>
    </li>
  </ol>
</details>


# 📝 À propos du projet
`aPhyloGeo` est un pipeline de bioinformatique dédié à l'analyse de la phylogéographie. `aPhyloGeo` est une application multiplateforme open source conçue par l'équipe du professeur Nadia Tahiri (Université de Sherbrooke, Québec, Canada). Elle est implémentée en Python. Cet outil peut être utilisé pour obtenir des arbres à partir de données climatiques des régions où les échantillons ont été collectés. Ces arbres climatiques sont ensuite utilisés pour une comparaison topologique et évolutive avec des arbres phylogénétiques issus d'alignements de séquences multiples (MSA) en utilisant la métrique des moindres carrés (LS) (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/). Les MSA qui produisent des arbres avec une valeur LS significative sont ensuite éventuellement enregistrés dans des dossiers avec leur arbre correspondant. Le fichier `output.csv` contient les informations de toutes les MSA significatives (voir la section Workflow pour plus de détails).

💡 Si vous utilisez notre algorithme dans votre recherche, veuillez citer notre article récent :
Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.
[Proceeding in SciPy 2022, Auxtin, TX, USA](https://conference.scipy.org/proceedings/scipy2022/pdfs/nadia_tahiri.pdf)

## Workflow

![](./img/Fig_1.png)

Figure 1: Le workflow de l'algorithme. Les opérations de ce workflow comprennent plusieurs blocs. Les blocs sont mis en évidence par trois couleurs différentes.
- **Le premier bloc** (couleur bleu clair) est responsable de la création des arbres basés sur les données climatiques - effectue la fonction de validation des paramètres d'entrée (voir le fichier YAML).
- **Le deuxième bloc** (couleur verte claire) est responsable de la création des arbres basés sur les données génétiques - effectue la fonction de validation des paramètres d'entrée (voir le fichier YAML).
- **Le troisième bloc** (couleur rose claire) permet la comparaison entre les arbres phylogénétiques (c'est-à-dire avec des données génétiques) et les arbres climatiques - étape de phylogéographie en utilisant la distance des moindres carrés (voir l'équation ci-dessus).

$$ls(T_1, T_2) = \sum_{1 \le i \le j \le n} \lvert \delta(i,j) - \xi(i,j) \rvert$$

où $T_1$ est l'arbre phylogénétique 1, $T_2$ est l'arbre phylogénétique 2, $i$ et $j$ sont deux espèces, $\delta(i,j)$ est la distance entre l'espèce $i$ et l'espèce $j$ dans $T_1$, $\xi(i,j)$ est la distance entre l'espèce $i$ et l'espèce $j$ dans $T_2$, et $n$ est le nombre total d'espèces.

Ceci est le bloc le plus important et la base de cette étude, à partir des résultats desquels l'utilisateur reçoit les données de sortie avec les calculs nécessaires.

De plus, notre approche est optimale car elle est élastique et s'adapte à n'importe quel ordinateur en utilisant le parallélisme et les GPU/CPUs disponibles en fonction de l'utilisation des ressources par unité de calcul (c'est-à-dire pour réaliser le traitement d'une seule fenêtre génétique - voir le flux de travail ci-dessous).
**Multiprocessing**: permet d'analyser plusieurs fenêtres simultanément (recommandé pour les grands ensembles de données)

Dans ce travail, nous avons utilisé des packages logiciels des versions suivantes: [Biopython](https://biopython.org/) version 1.79 (licence BSD 3 clauses), [Bio](https://pandas.pydata.org/) version 1.5.2 (licence New BSD) et [numpy](https://numpy.org/) version 1.21.6 (licence BSD 3 clauses).

# ⚒️ Installation

## Versions Linux UNIX et Mac OS
`aPhyloGeo` est disponible en tant que script Python.

### Prérequis
Ce package utilise l'outil de gestion des dépendances et d'emballage ```Poetry``` pour Python. Vous pouvez trouver le guide d'installation de Poetry ici : [Installation de Poetry](https://python-poetry.org/docs/#installation)

⚠️ Pour l'installation sous Windows, il est recommandé de lancer PowerShell en mode **Administrateur**.

Une fois Poetry installé, vous pouvez cloner le référentiel et installer le package à l'aide des commandes suivantes :

```
poetry install
```

### Utilisation
Poetry gérera automatiquement l'environnement virtuel. Si vous souhaitez utiliser l'environnement virtuel manuellement, vous pouvez utiliser la commande suivante :

```
poetry shell
```
⚠️ En supposant que Python 3.8 ou une version supérieure est installée sur la machine, le script devrait fonctionner correctement avec les bibliothèques installées.

Vous pouvez également lancer le package en utilisant la commande `make` depuis votre terminal lorsque vous êtes dans le `root`. Cette commande utilisera le fichier `Makefile` pour exécuter le script. Si vous utilisez la commande `make clean`, elle effacera le fichier `output.csv` précédemment créé avec la première commande.

# 🚀 Paramètres
Le logiciel `aPhyloGeo` peut être encapsulé dans d'autres applications et appliqué à d'autres données en fournissant un fichier YAML. Ce fichier inclura un ensemble de paramètres pour une manipulation facile.

- **Seuil de bootstrap**:  Nombre de réplicats seuil à générer pour chaque sous-MSA (chaque position de la fenêtre glissante)
- **Longueur de la fenêtre**: Taille de la fenêtre glissante
- **Pas**: Pas d'avancement de la fenêtre glissante
- **Choix de distance**: La distance des moindres carrés (LS) (version 1.0) sera étendue à la [métrique Robinson-Foulds (RF)](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)
- **Seuil de distance des moindres carrés**: Seuil de distance LS pour lequel les résultats sont les plus significatifs

# 📁 Example

## Description
Nous avons sélectionné seulement 5 des 38 lignées ayant des caractéristiques régionales pour une étude plus approfondie (voir Koshkarov et al., 2022). À partir des informations de localisation, des données de séquençage complet des nucléotides pour ces 5 lignées ont été collectées sur le [site NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). En cas de disponibilité de plusieurs résultats de séquençage pour la même lignée dans le même pays, nous avons sélectionné la séquence dont la date de collecte était la plus proche de la première date présentée. S'il y a plusieurs résultats de séquençage pour le même pays à la même date, la séquence avec le moins de caractères ambigus (N par nucléotide) est sélectionnée.

Bien que la sélection des échantillons se soit basée sur le cluster phylogénétique de la lignée et la transmission, la plupart des sites impliqués représentent des conditions météorologiques différentes. Comme le montre la Figure 2, les 5 échantillons impliquaient des températures allant de -4 C à 32,6 C, avec une température moyenne de 15,3 C. L'humidité spécifique variait de 2,9 g/kg à 19,2 g/kg avec une moyenne de 8,3 g/kg. La variabilité de la vitesse du vent et de l'irradiance descendante à courte longueur d'onde de la surface totale du ciel était relativement faible entre les échantillons par rapport à d'autres paramètres. La vitesse du vent variait de 0,7 m/s à 9,3 m/s avec une moyenne de 4,0 m/s, et l'irradiance descendante à courte longueur d'onde de la surface totale du ciel variait de 0,8 kW-hr/m2/jour à 8,6 kW-hr/m2/jour avec une moyenne de 4,5 kW-hr/m2/jour. Contrairement aux autres paramètres, 75% des villes impliquées reçoivent moins de 2,2 mm de précipitations par jour, et seules 5 villes ont plus de 5 mm de précipitations par jour. La précipitation minimale est de 0 mm/jour, la précipitation maximale est de 12 mm/jour, et la valeur moyenne est de 2,1 mm/day.

## Entrée

L'algorithme prend en entrée deux fichiers avec les définitions suivantes :

- 🧬 **Genetic file** avec l'extension fasta. Le premier fichier ou ensemble de fichiers contiendra les informations de séquence génétique des ensembles d'espèces sélectionnés pour l'étude. Le nom du fichier doit permettre de connaître le nom du gène. Il est donc fortement recommandé de suivre la nomenclature suivante : gene_name.fasta.
- ⛅ **Climatic file** avec l'extension csv. Le deuxième fichier contiendra les informations d'habitat pour les ensembles d'espèces sélectionnés pour l'étude. Chaque ligne représentera l'identificateur de l'espèce et chaque colonne représentera une condition climatique.

## Sortie
L'algorithme retournera un fichier csv contenant des informations de toutes les MSAs pertinentes (voir la section Worflow pour plus de détails). Les fenêtres glissantes d'intérêt sont celles avec un support bootstrap intéressant (c'est-à-dire indiquant la robustesse de l'arbre) et une similarité élevée à la condition climatique en question (c'est-à-dire basée sur la valeur `LS`). Ils indiqueront, entre autres, le nom du gène, la position du début et de la fin de la fenêtre glissante, la valeur moyenne du bootstrap, la valeur LS et enfin la condition climatique pour laquelle cette zone génétique expliquerait l'adaptation des espèces dans un environnement donné.


# ✔️ Références

1️⃣ Calcul de la distance entre l'arbre phylogénétique: `Métrique des moindres carrés`
+ [Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis. Models and estimation procedures. American journal of human genetics, 19(3 Pt 1), 233.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/)
+ [Felsenstein, J. (1997). An alternating least squares approach to inferring phylogenies from pairwise distances. Systematic biology, 46(1), 101-111.](https://pubmed.ncbi.nlm.nih.gov/11975348/)
+ [Makarenkov, V., & Lapointe, F. J. (2004). A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. Bioinformatics, 20(13), 2113-2121.](https://pubmed.ncbi.nlm.nih.gov/15059836/)

2️⃣ Calcul de la distance entre l'arbre phylogénétique: `Métrique Robinson-Foulds`
+ [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)

3️⃣ Description complète de l'ensemble de données : `Analyse des données génétiques et climatiques de SARS-CoV-2`
+ [Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.](https://conference.scipy.org/proceedings/scipy2022/nadia_tahiri.html)

# 📧 Contact
Veuillez nous envoyer un courriel à l'adresse suivante: <Nadia.Tahiri@USherbrooke.ca> pour toute question ou commentaire.
