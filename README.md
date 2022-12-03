﻿﻿﻿﻿﻿﻿﻿<h1  align="center"> aPhyloGeo <p align='center'> 
        [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
        [![Contributions](https://img.shields.io/badge/contributions-welcome-blue.svg)](https://pysd.readthedocs.io/en/latest/development/development_index.html)
        [![Py version](https://img.shields.io/pypi/pyversions/pysd.svg)](https://pypi.python.org/pypi/pysd/)
        [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftahiri-lab%2FaPhylogeo&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
        [![GitHub release](https://img.shields.io/github/v/release/tahiri-lab/aPhylogeo.svg?maxAge=3600)](https://github.com/tahiri-lab/aPhylogeo/releases/)
        </p>


<h2  align="center">Multi-platform application for analyze phylogenetic trees with climatic parameters</h2>

<details open>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About the project</a>
    </li>
    <li>
      <a href="#Installation">Installation</a>
      <ul>
        <li><a href="#Linux-UNIX-and-Mac-OS-versions">Linux/UNIX and Mac OS versions</a></li>
      </ul>
    </li>
    <!--<li> Available analyses</li>
      <ul>
        <li><a href="#Group-creation">Group creation</a></li>
        <li><a href="#SimPlot-analysis">SimPlot analysis</a></li>
        <li><a href="#Similarity-networks">Similarity networks</a></li>
        <li><a href="#BootScan-analysis">BootScan analysis</a></li>
        <li><a href="#Findsites">Findsites</a></li>
        <li><a href="#Detection-of-recombination">Detection of recombination</a></li>
      </ul>-->
    <li>
      <a href="#Reference">Reference</a>
    </li>
    <li>
      <a href="#contact">Contact</a>
    </li>
  </ol>
</details>


# About the project

Link to the French version: [French version](https://github.com/tahiri-lab/aPhylogeo/blob/main/README_fr.md)

`aPhylogeo` is a bioinformatics pipeline dedicated to the analysis of phylogeography. This tool can be used to obtain trees from climatic data of the regions where the samples have been collected. Those climatic trees are then used for topological comparison against phylogenetic trees from multiple sequence alignments (MSAs) using the [Robinson-Foulds (RF) metric](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub). MSAs that yield trees with a significant `RF` value are then saved in folders with their respective tree. The `output.csv` file contains the informations of all the significant MSAs informations.

## Workflow

![](./img/workflow_en.png)


The workflow of the algorithm. The operations within this workflow include several blocks. The blocks are highlighted by three different colors. The first block (the light pink color) is responsible for creating the trees based on the climate data. The second block (the dark yellow color) performs the function of input parameter validation. The third block (the light-yellow color) allows the creation of phylogenetic trees. This is the most important block and the basis of this study, through the results of which the user receives the output data with the necessary calculations.


# Installation

## Linux UNIX and Mac OS versions
SimPlot++ is available as a Python script.

### Prerequisites
Before using this program, make sure that you have installed all the necessary libraries for it to work properly. To do this, simply type the following command:

```
pip install -r requirements.txt
```

### Python script
A `requirements.txt` file containing all required libraries is available in the GitHub repository.

Assuming Python 3.8 or higher is installed on the machine, the script should run well with the libraries installed.

<u>Here is an example of how to run the script in Linux/UNIX or Mac OS:</u>
1. After downloading the source code, go to the folder containing `main.py`.
2. If you do not have `virtualenv` installed, run `python3 -m pip install --user virtualenv`
3. Create a new virtual environment (venv) in your terminal using `python3 -m venv aPhyloGeo_env`.
4. Still in the terminal, enter the new venv using `source aPhyloGeo_env/bin/activate`.
5. Install the required libraries using `python3 -m pip install -r requirements.txt`.
6. Launch SimPlot++ using `python3 main.py`.



### Creation of climatic trees

In order to obtain the climatic data trees, you need to create a csv file containing the data to study ( see the example [data.csv](./data/data.csv) ).

This file should have a structure similar to the example file. A column containing the names of the specimens should be present since the name of this column will be requested when the program is launched. **The column names should not contain parentheses**!

The other columns should contain the climatic variables that are to be studied. The values must be **numeric** for the program to work properly.

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

Make sure that the column names are identical to the ones you entered so that the program can find them. The names of the sequences to be studied and the names of the specimens in this file must also match exactly.

The climatic trees should then be created in the current folder with the column names followed by ``newick`` (e.g. ``T_max_to_2m_C_newick``). 


### Creating phylogenetic trees and phylogeographic analysis

1. The `.fasta` sequence files must be put in the folder [data](./data). Note that the header of the sequence file must respect this format (Note: The sequences obtained on [GISAID](https://www.gisaid.org) all respect this format):

```>hCoV-19/Nigeria/S38/2020|EPI_ISL_2399462|2020-10-12```

In order to be able to differentiate each sequence with a different name, this tool will isolate the 3rd word encountered in the header. In this example, the sequence will be named `S38`. In case this word is longer than 10 characters, the program will only take the first ten : 

```>hCoV-19/India/GJ-GBRC560b/2021|EPI_ISL_1677798|2021-01-12```

In this case, the sequence name will only be `GJ-GBRC560` since it is exactly 10 characters long. This is important to understand, especially for the names you will give to the sequences in the `csv` climatic data file, which must be identical to the single names.

2. Once the sequences have been placed in the [data](./data) folder, all you have to do is run the program with `make` command. Here is an example of what should appear on the terminal and an example of the entries:

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
Depending on the bootstrap value and the `Robinson-Fould metric` found for each MSA ( Multiple Sequence Alignment ), the climatic tree of interest with their corresponding MSA will be kept in the folders of each selected [gene](./output) or in the [reference](./output/reference_gene) folder if the complete sequence has been studied. Finally, the [output.csv](output.csv) file will contain all the significant MSAs and their related informations such as their position on the sequence, the gene related to this position, their bootstrap value and their `RF metric`.

In this work, we applied software packages of the following versions: MUSCLE version 3.18 (GNU GENERAL PUBLIC LICENSE), PHYLIP version 3.18 (open source license), RAxML version 8.2.12 (GNU GENERAL PUBLIC LICENSE).

## Potential problems encountered

+ For `macOS` users, it is likely that your computer is blocking access to the `MUSCLE` program. If this is the case, go to the privacy settings on your machine and give the program access

## References

1. Sequence alignment tool : `MUSCLE`
+ [Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797.](https://academic.oup.com/nar/article/32/5/1792/2380623)
doi:10.1093/nar/gkh340]
+ [Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity BMC Bioinformatics, (5) 113.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113)
doi:10.1186/1471-2105-5-113

2. Phylip Package: `Seqboot`, `DNADist`, `Neighbor`, `Consense`
    [J. Felsenstein. 1989. PHYLIP - Phylogeny Inference Package (Version 3.2) . Cladistics. 5: 164-166.](https://evolution.genetics.washington.edu/phylip.html)

3. Calculation of distance between phylogenetic tree: `Robinson-Foulds metric`
    [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)

4. Phylogenetic analysis: `RAxML`
    [A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053?login=true)


# Contact
Please email us at : <Nadia.Tahiri@USherbrooke.ca> for any question or feedback.
