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

# Abstract
*aPhyloGeo* is a Python library that explores the relationship between species evolution and environmental pressures, particularly climate. By integrating genetic and climatic data, it helps researchers investigate evolutionary adaptations and identify genetic regions influenced by environmental factors.

The core feature of the software is a comprehensive phylogenetic analysis pipeline with three investigation levels: 1) genetic relationships, 2) climatic impact assessment, and 3) biogeographic correlations. This approach facilitates understanding how species adapt to their environments. 

*aPhyloGeo* employs algorithms using metrics like least squares [@felsenstein1997alternating], Euclidean, and Robinson-Foulds [@robinson1981comparison] distances to ensure statistically sound correlations. Its modular structure offers flexibility, and its open-source nature promotes collaboration.

Available as a PyPI package, *aPhyloGeo* enhances understanding of evolutionary processes and informs conservation efforts, helping prioritize species and populations for preservation in the face of climate change.

# Statement of Need

The rapid impacts of climate change and anthropogenic variables on biodiversity and population dynamics underscore the necessity for more advanced tools capable of resolving the complexities of ecosystems under perturbation. Biologists utilize phylogeographic approaches to closely examine the interplay between the genetic structures of study populations and their geographic distributions, considering both current and historical geoclimatic contexts.

This software package, *aPhyloGeo*, is designed for detailed phylogeographic analysis (i.e., Genetic dataset and climatic dataset).
Given the urgency of the current climate crisis and the anticipated future challenges, there is a pressing need to develop tools that not only meet but also exceed bioinformatics software development standards. These tools will be crafted to enable accurate characterization of genetic diversity and phenotypic traits in strict accordance with environmental conditions, thus empowering researchers to address critical questions in the field, such as:

- How have past climatic fluctuations shaped the current patterns of genetic diversity and geographic distribution within species?
- How will species distribution evolve under different future climate scenarios? Is it possible to identify potential refuges where species could persist?
- Which populations are most vulnerable to climate-related extinction based on their genetic diversity and adaptive potential?
- Are there genetic signatures of local adaptation that reveal how populations have evolved in response to specific environmental pressures?

This research aims to make a significant contribution to our understanding of the evolving ecological landscape and provide the scientific community with robust tools for comprehensive analysis and interpretation. *aPhyloGeo* will enable researchers to unravel the complex interplay between genetics, geography, and environment, informing conservation strategies and predicting the impacts of ongoing environmental change.

# State of the Field - Advancements in Genomic Analysis
Advances in genomic analysis, notably in identifying gene sub-sequences [@nadia_tahiri-proc-scipy-2022], and applying them to SARS-CoV-2 data [@nadia_tahiri-proc-scipy-2023], enhance understanding of genetic adaptation across species and environments. In phylogeography, though various Python packages provide easy analysis, none seamlessly integrates genetic and environmental data. While Biopython handles genetic sequences well [@cornish2021biopython], it lacks environmental integration. DendroPy visualizes phylogeographic patterns but needs additional tools (actually our team works on [iPhyloGeo++](https://github.com/tahiri-lab/iPhyloGeo_plus_plus), a software). SciPy offers statistical utilities but requires expertise. GeoPandas manages geospatial data but lacks genetic integration. A user-friendly Python package specifically for phylogeography is needed.

Statistical methods like generalized linear models and mixed models help quantify the relationship between genetic variation and environmental factors. The evolution of these tools alongside high-throughput sequencing and environmental data availability creates opportunities in evolutionary biology, ecology, and conservation. *aPhyloGeo* integrates genetic and climatic data, aiming to enhance understanding of biodiversity dynamics.

# Pipeline

Navigating the *aPhyloGeo* workflow (refer to \autoref{fig:figure1}) is indispensable to fully harness the potential of this bioinformatics pipeline. The visual representation in \autoref{fig:figure1} outlines the key steps for conducting phylogeographic analysis with optimized effectiveness.
For a more detailed understanding, please refer to the comprehensive tutorial [provided with the software](https://github.com/tahiri-lab/aPhyloGeo/wiki/Tutorial).

![The workflow of the algorithm. The operations within this workflow include several blocks.\label{fig:figure1}](../img/workflow_en.png)

The diagram below illustrates the workflow of the algorithm, consisting of several key blocks, each highlighted with a distinct color (refer to the [wiki page](https://github.com/tahiri-lab/aPhyloGeo/wiki/Worflow).

## Multiprocessing

The algorithm supports multiprocessing, allowing the simultaneous analysis of multiple sliding windows (i.e., corresponds to the alignment sub-sequence where the starting point is located after the first position of the alignment, and the ending point is located before the last position of the alignment) within the genetic data. This feature is particularly recommended for large datasets, as it significantly speeds up the analysis by dividing the input sequences into smaller chunks that can be processed in parallel.

## Dependencies

This work relies on the following main software packages:

- [ete3](https://pypi.org/project/ete3/) version 3.1.3 [@huerta2016ete].
- [Bio](https://pypi.org/project/bio/) version 1.5.9 [@cock2009biopython].
- [robinson-foulds](https://pypi.org/project/robinson-foulds/) version 1.2 [@huerta2016ete].
- [dendropy](https://pypi.org/project/DendroPy/) version 4.6.1 [@sukumaran2010dendropy].

# Methods

## Tree Comparison

In the comparison of phylogenetic trees, which are constructed based on genetic data, with climatic trees, a crucial step involves applying a phylogeography approach. This includes the utilization of Robinson and Foulds distance for topology evaluation and Least Squares distance for assessing branch length differences.

## Editing Multiple Sequence Alignment Methods

Multiple Sequence Alignment (MSA) holds immense significance in bioinformatics as it serves as a foundational step for the comparison and analysis of biological sequences. Here is an in-depth overview of some widely used MSA methods: 1) **Pairwise Alignment** [@li2018minimap2], 2) **MUSCLE** [@edgar2004muscle], 3) **CLUSTALW** [@hung2016sequence], and **MAFFT** [@katoh2013mafft].

## Similarity Methods

Sequences with notable variability were specifically retained for analysis. The dissimilarity assessment between each sequence pair involved the application of an extensive set of 8 metrics: 1) **Hamming distance** [@labib2019hamming], 2) **Levenshtein distance** [@yujian2007normalized], 3) **Damerau-Levenshtein distance** [@zhao2019string], 4) **Jaro similarity** [@pradhan2015review], 5) **Jaro-Winkler similarity** [@pradhan2015review], 6) **Smith–Waterman similarity** [@waterman1978similarity], 7) **Jaccard similarity** [@bag2019efficient], and 8) **Sørensen-Dice similarity** [@li2020generic].

This comprehensive methodology ensures a nuanced and high-quality analysis, contributing to a deeper understanding of sequence distinctions.


# Conclusion

The *aPhyloGeo* pipeline serves as an integrative framework, bringing together a variety of advanced analytical methodologies for diverse datasets, covering both genetic and climatic aspects. By consolidating these analyses within a unified platform, users can simplify their exploration of different tools while ensuring greater reproducibility in research outcomes.

Looking ahead, *aPhyloGeo* aims to integrate new functionalities, including clustering techniques based on similarity derived from multiple sequence alignments and a more computationally efficient alignment methodology. The incorporation of novel metrics, such as the Euclidean distance and Robinson-Foulds distance, aims to provide users with improved insight regarding their datasets through a comprehensive assessment of genetic diversity.

By adhering to best practices in software development and embracing open-source principles, *aPhyloGeo* not only provides a reliable and adaptable platform for current phylogeographic research but also lays the groundwork for future expansion and innovation. The ongoing integration of cutting-edge methodologies, such as clustering based on multiple sequence alignments and optimized alignment algorithms, demonstrates our commitment to continually enhance the capabilities of aPhyloGeo. Moreover, the incorporation of novel metrics like the Euclidean distance and Robinson-Foulds distance empowers users with enhanced decision-making tools for assessing genetic diversity and its relationship to environmental factors. Through these targeted improvements, aPhyloGeo aims to foster reproducible, accessible, and comprehensive phylogeographic analysis, ultimately deepening our understanding of the complex interactions between species and their environments.

# Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada, Fonds de recherche du Québec – Nature et technologies, the University of Sherbrooke grant, and the Centre de recherche en écologie de l'UdeS (CREUS). The computations were performed on resources provided by Compute Canada and Compute Quebec, the national and provincial infrastructure for high-performance computing and data storage. The authors would like to thank the students of the University of Sherbrooke and the Université du Québec à Montréal for their significant contributions to the development of the software. Finally, the authors would like to thank the reviewers and the editor for their valuable comments.

# References
