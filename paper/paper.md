---
title: 'aPhyloGeo: a multi-platform Python package for analyzing phylogenetic trees with climatic parameters'
tags:
  - bioinformatics
  - consensus
  - metrics
  - multiple sequence alignment
  - phylogeny
  - phylogeography
authors:
  - name: My-Linh Luu
    affiliation: 2
    corresponding: false
    email: mylinh_93@hotmail.com
  - name: Georges Marceau
    affiliation: 2
    orcid: 0009-0005-9827-5735
    corresponding: false
    email: gmaccounts@icloud.com
  - name: David Beauchemin
    affiliation: 2
    orcid: 0009-0009-1338-9012
    corresponding: false
    email: david.beauchemin.3@courrier.uqam.ca
  - name: Nadia Tahiri
    affiliation: 1
    orcid: 0000-0002-1818-208X
    corresponding: true
    email: Nadia.Tahiri@USherbrooke.ca
affiliations:
  - name: département d’Informatique, Université de Sherbrooke, 2500 Boulevard de l’Université, Sherbrooke, Québec J1K 2R1, Canada
    index: 1
  - name: département d’Informatique, Université du Québec à Montréal, 201, avenue du Président-Kennedy, Montréal, Québec H2X 3Y7, Canada
    index: 2
date: 02 March 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
*aPhyloGeo*, a versatile and open-source Python application available at [PyPI](https://pypi.org/project/aphylogeo/), is designed to elucidate the complex relationship between species evolution and environmental pressures, with a particular focus on climate. By integrating genetic and climatic data, *aPhyloGeo* empowers researchers to investigate the mechanisms of evolutionary adaptation and pinpoint genetic regions potentially influenced by environmental factors.

The software's core strength lies in its comprehensive phylogenetic analysis pipeline, encompassing three distinct levels of investigation: genetic relationships, climatic impact assessment, and biogeographic correlations. This multi-faceted approach facilitates a holistic understanding of how species evolve and adapt to their environments. For example, researchers studying the genetic basis of high-altitude adaptation in birds could utilize aPhyloGeo to construct phylogenetic trees from genetic data, analyze oxygen levels across different altitudes, and identify correlations between specific genes and hypoxic conditions.

In another scenario, scientists investigating the impact of climate change on marine biodiversity could employ *aPhyloGeo* to examine the genetic diversity of coral species, assess changes in sea surface temperatures over time, and pinpoint genetic markers associated with thermal tolerance.
These examples demonstrate the wide range of research questions that *aPhyloGeo* can address, making it an invaluable tool for evolutionary biologists, ecologists, and conservationists alike.

Underlying *aPhyloGeo*'s analyses are robust algorithms employing metrics such as least squares distance [@felsenstein1997alternating] , Euclidean distance, and Robinson-Foulds distance [@robinson1981comparison] to quantify similarity across different levels. This rigorous approach ensures that the identification of correlations is statistically sound, while adhering to the principles of phylogenetic inference [@gascuel2006neighbor]. The software's modular design and Python interface offer flexibility, allowing users to tailor analyses to their specific research questions and datasets. Additionally, *aPhyloGeo*'s open-source nature fosters collaboration and transparency within the scientific community.

By enabling researchers to explore the complex interplay between genetics and environment, *aPhyloGeo* contributes to a deeper understanding of evolutionary processes. This knowledge not only enhances our appreciation of the natural world but also informs conservation efforts in the face of climate change and other environmental challenges. By identifying genetic adaptations to changing environments, *aPhyloGeo* can help prioritize species and populations for conservation, ultimately contributing to the preservation of biodiversity on our planet.

# Statement of Need

The rapid impacts of climate change and anthropogenic variables on biodiversity and population dynamics underscore the necessity for more advanced tools capable of resolving the complexities of ecosystems under perturbation. Biologists utilize phylogeographic approaches to closely examine the interplay between the genetic structures of study populations and their geographic distributions, considering both current and historical geoclimatic contexts.

This software package, aPhyloGeo, is dedicated to advancing state-of-the-art bioinformatics tools specifically designed for detailed phylogeographic analysis. Given the urgency of the current climate crisis and the anticipated future challenges, there is a pressing need to develop tools that not only meet but also exceed bioinformatics software development standards. These tools will be crafted to enable accurate characterization of genetic diversity and phenotypic traits in strict accordance with environmental conditions, thus empowering researchers to address critical questions in the field, such as:
- How have past climatic fluctuations shaped the current patterns of genetic diversity and geographic distribution within species?
- How will species distribution evolve under different future climate scenarios? Is it possible to identify potential refuges where species could persist?
- Which populations are most vulnerable to climate-related extinction based on their genetic diversity and adaptive potential?
- Are there genetic signatures of local adaptation that reveal how populations have evolved in response to specific environmental pressures?

By maintaining the highest standards, this research aims to make a significant contribution to our understanding of the evolving ecological landscape and provide the scientific community with robust tools for comprehensive analysis and interpretation. aPhyloGeo will enable researchers to unravel the complex interplay between genetics, geography, and environment, informing conservation strategies and predicting the impacts of ongoing environmental change.

# State of the Field - Advancements in Genomic Analysis
The field of genomic analysis has progressed significantly in recent years, notably in the creation of tools and algorithms to explore the intricate relationship between genetic variation and environmental factors. The Tahiri lab team's innovative 2021 algorithm for identifying sub-sequences within genes [@nadia_tahiri-proc-scipy-2022], and its subsequent application to SARS-CoV-2 data in 2023 [@nadia_tahiri-proc-scipy-2023], stand out as substantial contributions to this field, enhancing our comprehension of the genetic underpinnings of adaptation across various species and environments.

In the broader field of phylogeography, substantial methodological advancements have also occurred. Several Python packages provide functionalities pertinent to phylogeographic analysis, but often in a fragmented way. Biopython [@cornish2021biopython], a cornerstone in bioinformatics, excels at handling genetic sequences and basic phylogenetic tasks, yet falls short in integrating environmental data. DendroPy, a robust library for phylogenetic trees, aids in visualizing phylogeographic patterns but requires additional tools for comprehensive analysis. While SciPy's statistical prowess could be harnessed for custom analyses, its complexity demands a strong background in statistical programming. GeoPandas, adept at handling geospatial data, is useful for mapping genetic or environmental distributions, but lacks seamless integration with genetic data analysis tools. In summary, while powerful individual tools exist, a comprehensive and user-friendly Python package specifically designed for phylogeographic analysis remains a gap to be filled.

Statistical approaches, including generalized linear models (GLMs) and mixed models, are increasingly used to investigate the relationship between genetic variation and environmental variables. These methods enable researchers to quantify the relative influence of various factors, such as climate, geography, and demography, on observed patterns of genetic diversity.

The continuous refinement of these tools and methodologies, coupled with the growing availability of high-throughput sequencing technologies and environmental data, has opened up exciting new research avenues in evolutionary biology, ecology, and conservation. aPhyloGeo builds upon these advancements, providing a unified platform for integrating genetic and climatic data to address a wide array of phylogeographic questions. By bridging the gap between genomics and environmental science, aPhyloGeo aims to contribute to a more comprehensive understanding of the forces shaping biodiversity in a changing world.

# Pipeline

Navigating the *aPhyloGeo* workflow (refer to \autoref{fig:figure1}) is indispensable to fully harness the potential of this bioinformatics pipeline. The visual representation in \autoref{fig:figure1} outlines the key steps for conducting phylogeographic analysis with optimized effectiveness.

![The workflow of the algorithm. The operations within this workflow include several blocks.\label{fig:figure1}](../img/workflow_en.png)

The diagram below illustrates the workflow of the algorithm, consisting of several key blocks, each highlighted with a distinct color.

- **First Block (Light Blue):** This block creates climate trees based on input climate data (CSV file) and validates the input parameters using a YAML file. More precisely, the climate trees were generated by calculating the pairwise differences between each value of the species' habitats, normalized between the minimum and maximum of the parameter. This process resulted in a symmetric square matrix. From this matrix, the climate tree was inferred using the Neighbor-Joining method. This involves processing climatic variables such as temperature, precipitation, and elevation to construct phylogenetic trees that represent the relationships between geographic locations based on their climatic similarity.

- **Second Block (Light Green):** This block creates phylogenetic trees based on input genetic data and performs input parameter validation (refer to the [YAML file](../aphylogeo/params.yaml)). This entails aligning DNA or amino acid sequences, inferring phylogenetic relationships using various methods (e.g., maximum likelihood, Bayesian inference), and assessing the statistical support for the inferred tree topology.

- **Third Block (Light Pink):**  The third block, referred to as the phylogeography step, is the crux of the analysis. It compares the genetic trees (representing evolutionary relationships) with the climate trees (representing environmental similarity). This comparison utilizes either the Robinson-Foulds distance or the Least Squares distance to quantify the degree of congruence between the two tree types.  The output of this step includes:
  - Topological congruence statistics: Quantifying the degree of similarity between the genetic and climate trees.
  - Co-phylogenetic visualizations: Graphical representations highlighting the associations between genetic lineages and climatic niches.
  - Statistical tests: Assessing the significance of the observed phylogeographic patterns.
This third block is pivotal, forming the basis from which users obtain output data (i.e., name of gene, name of climate parameter, bootstrap value, Robinson-Foulds distance, entropy distance, least-square distance, the starting position and the ending position of windows, and climatic and genetic trees) with essential calculations (i.e., distances, tree inference, sequence alignment). Our approach is optimized to adapt to various computing environments through elasticity and utilize parallelism and available GPUs/CPUs based on resource usage per unit of computation. This flexibility enables efficient processing of a single genetic window, as outlined in the workflow below.

## Multiprocessing

The algorithm supports multiprocessing, allowing the simultaneous analysis of multiple sliding windows (i.e., corresponds to the alignment sub-sequence where the starting point is located after the first position of the alignment, and the ending point is located before the last position of the alignment) within the genetic data. This feature is particularly recommended for large datasets, as it significantly speeds up the analysis by dividing the input sequences into smaller chunks that can be processed in parallel.

## Dependencies

This work relies on the following main software packages:

- [ete3](https://pypi.org/project/ete3/) version 3.1.3 [@huerta2016ete], available under the GNU General Public License (GPL) (GPLv3), is used for phylogenetic tree manipulation and visualization.
- [Bio](https://pypi.org/project/bio/) version 1.5.9 [@cock2009biopython], available under the New BSD License, provides a wide range of bioinformatics tools and functionalities.
- [robinson-foulds](https://pypi.org/project/robinson-foulds/) version 1.2 [@huerta2016ete], available under the GNU General Public License v3 (GPLv3), is utilized for calculating the Robinson-Foulds distance between phylogenetic trees.
- [dendropy](https://pypi.org/project/DendroPy/) version 4.6.1 [@sukumaran2010dendropy], available under the BSD License (BSD), is employed for phylogenetic tree manipulation and analysis.
# Methods

## Tree Comparison

In the comparison of phylogenetic trees, which are constructed based on genetic data, with climatic trees, a crucial step involves applying a phylogeography approach. This includes the utilization of Robinson and Foulds distance for topology evaluation and Least Squares distance for assessing branch length differences.

## Editing Multiple Sequence Alignment Methods

Multiple Sequence Alignment (MSA) holds immense significance in bioinformatics as it serves as a foundational step for the comparison and analysis of biological sequences. Here is an in-depth overview of some widely used MSA methods:

- **Pairwise Alignment**: Fundamental in comparing two sequences [@li2018minimap2].
- **MUSCLE**: Multiple Sequence Comparison by Log-Expectation, a popular tool for high-quality MSA [@edgar2004muscle].
- **CLUSTALW**: A widely-used software for multiple sequence alignment [@hung2016sequence].
- **MAFFT**: Multiple Alignment using Fast Fourier Transform, known for its accuracy and efficiency [@katoh2013mafft].

## Similarity Methods

To enhance the algorithm's performance, a meticulous approach was adopted. Sequences with notable variability were specifically retained for analysis. The dissimilarity assessment between each sequence pair involved the application of an extensive set of 8 metrics:

1. **Hamming distance**: Measures the difference between two strings of equal length [@labib2019hamming].
2. **Levenshtein distance**: Evaluates the minimum number of single-character edits required to transform one sequence into another [@yujian2007normalized].
3. **Damerau-Levenshtein distance**: Similar to Levenshtein distance, with an additional operation allowing transpositions of adjacent characters [@zhao2019string].
4. **Jaro similarity**: Computes the similarity between two strings, considering the number of matching characters and transpositions [@pradhan2015review].
5. **Jaro-Winkler similarity**: An enhancement of Jaro similarity, giving more weight to common prefixes [@pradhan2015review].
6. **Smith–Waterman similarity**: Utilizes local sequence alignment to identify similar regions within sequences [@waterman1978similarity].
7. **Jaccard similarity**: Measures the similarity between finite sample sets [@bag2019efficient].
8. **Sørensen-Dice similarity**: Particularly useful for comparing the similarity of two samples [@li2020generic].

This comprehensive methodology ensures a nuanced and high-quality analysis, contributing to a deeper understanding of sequence distinctions.


# Conclusion

The *aPhyloGeo* pipeline serves as an integrative framework, bringing together a variety of advanced analytical methodologies for diverse datasets, covering both genetic and climatic aspects. By consolidating these analyses within a unified platform, users can simplify their exploration of different tools while ensuring greater reproducibility in research outcomes.

Looking ahead, *aPhyloGeo* aims to integrate new functionalities, including clustering techniques based on similarity derived from multiple sequence alignments and a more computationally efficient alignment methodology. The incorporation of novel metrics, such as the Euclidean distance and Robinson-Foulds distance, aims to provide users with improved insight for making nuanced decisions regarding their datasets through a comprehensive assessment of genetic diversity.

By adhering to best practices in software development and embracing open-source principles, *aPhyloGeo* not only provides a reliable and adaptable platform for current phylogeographic research but also lays the groundwork for future expansion and innovation. The ongoing integration of cutting-edge methodologies, such as clustering based on multiple sequence alignments and optimized alignment algorithms, demonstrates our commitment to continually enhance the capabilities of aPhyloGeo. Moreover, the incorporation of novel metrics like the Euclidean distance and Robinson-Foulds distance empowers users with enhanced decision-making tools for assessing genetic diversity and its relationship to environmental factors. Through these targeted improvements, aPhyloGeo aims to foster a new era of reproducible, accessible, and comprehensive phylogeographic analysis, ultimately deepening our understanding of the complex interactions between species and their environments.

# Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada, Fonds de recherche du Québec - Nature et technologie, the University of Sherbrooke grant, and the Centre de recherche en écologie de l'UdeS (CREUS). The author would like to thank the Department of Computer Science, University of Sherbrooke, Quebec, Canada for providing the necessary resources to conduct this research. The computations were performed on resources provided by Compute Canada and Compute Quebec - the National and Provincial Infrastructure for High-Performance Computing and Data Storage. The author would like to thank the students of the University of Sherbrooke and the Université du Québec à Montréal for their great contribution to the development of the software.

# References
