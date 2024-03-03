---
title: 'aPhyloGeo: Multi-platform application for analyze phylogenetic trees with climatic parameters'
tags:
  - bioinformatics
  - consensus
  - metrics
  - multiple sequence alignment
  - phylogeny
  - phylogeography
authors:
  - name: Nadia Tahiri
    affiliation: 1
    orcid: 0000-0002-1818-208X
    corresponding: true
    email: Nadia.Tahiri@USherbrooke.ca
affiliations:
  - name: département d’Informatique, Université de Sherbrooke, 2500 Boulevard de l’Université, Sherbrooke, Québec J1K 2R1, Canada
    index: 1
date: 02 March 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
The cross-platform application for phylogenetic tree analysis with climate parameters, `aPhyloGeo`, is a robust pipeline designed for comprehensive phylogenetic analyses using genetic and climate data. This Python API, available on [PyPI](https://pypi.org/project/aphylogeo/), offers a suite of analyses tailored to various scenarios, enabling the examination of datasets at three distinct levels: 1) genetic, 2) climatic, and 3) biogeography correlation, all within a unified package. Similarity at these levels, evaluated through metrics such as least squares distance and Robinson and Foulds distance, significantly influences the assumptions guiding the identification of correlations between a species' genetics and its habitat during the reconstruction of the multiple alignment necessary for phylogenetic inference.

By utilizing the `aPhyloGeo` Python API, users can programmatically implement sophisticated phylogenetic analyses without the need for a graphical interface. This API provides a powerful and flexible toolset for conducting analyses, allowing users to tailor the application to their specific research needs. Through this approach, aPhyloGeo facilitates a nuanced understanding of the interplay between genetic evolution and environmental factors in the context of species adaptation, all within the Python programming environment.

By selecting an appropriate gene list for the available data defined on a set of species to explain the adaptation of the species according to the Darwinian hypothesis, the user can be confident that these assumptions are taken into account in `aPhyloGeo`.

# Statement of Need

The rapid impacts of climate change and anthropogenic variables on biodiversity and population dynamics underscore the need for more advanced tools capable of resolving the complexities of ecosystems under perturbation. Biologists use phylogeographic approaches to closely examine the interplay between the genetic structures of study populations and their geographic distributions, taking into account both current and historical geoclimatic contexts.

This software package is dedicated to advancing state-of-the-art bioinformatics tools specifically designed for detailed phylogeographic analysis. Given the urgency of the current climate crisis (COP27 - Climate Change and COP15 - Convention on Biological Diversity) and the anticipated future challenges, there is an urgent need to develop tools that not only meet but exceed bioinformatics software development standards. These tools will be designed to allow accurate characterization of genetic diversity and phenotypic traits in strict accordance with environmental conditions. By maintaining the highest standards, this research aims to make a significant contribution to our understanding of the evolving ecological landscape and provide the scientific community with robust tools for comprehensive analysis and interpretation.

# State of the Field

In 2021, Tahiri lab team [@nadia_tahiri-proc-scipy-2022] proposed a new algorithm to allow finding sub-sequences of genes giving an increased topological similarity between the reference tree (obtained from gene sequences) and the phylogenetic tree (obtained from genome sequences). It can help find which genes or subparts of a gene are sensitive or favourable to a given environment. Finally, [@nadia_tahiri-proc-scipy-2023] explored this algorithm with SARS-CoV-2 data.

# Pipeline

Exploring the aPhyloGeo workflow (\autoref{fig:figure1}) is essential to harness the full potential of this bioinformatics pipeline. Follow these steps to perform phylogeographic analysis effectively:

## Algorithm Workflow

![The workflow of the algorithm. The operations within this workflow include several blocks.\label{fig:figure1}](../img/workflow_en.png)

The diagram below illustrates the workflow of the algorithm, consisting of several key blocks, each highlighted with a distinct color.

- **First Block (Light Blue):** This block is responsible for creating trees based on climate data and performs input parameter validation (refer to the YAML file).

- **Second Block (Light Green):** This block focuses on creating trees based on genetic data and conducts input parameter validation (refer to the YAML file).

- **Third Block (Light Pink):** The third block facilitates the comparison between phylogenetic trees (genetic data) and climatic trees, denoted as the phylogeography step. It utilizes the Robison and Foulds distance or Least Square distance.

This third block is pivotal to the study, forming the basis from which users obtain output data with essential calculations. Our approach is optimal, adapting to various computing environments through elasticity and utilizing parallelism and available GPUs/CPUs based on resource usage per unit of computation. This flexibility enables efficient processing of a single genetic window, as outlined in the workflow below.

## Multiprocessing

The algorithm supports multiprocessing, allowing simultaneous analysis of multiple windows. This feature is particularly recommended for large datasets.

## Dependencies

This work relies on the following software packages:

- [Biopython](https://biopython.org/) version 1.79 (BSD 3-Clause License)
- [Bio](https://pandas.pydata.org/) version 1.5.2 (New BSD License)
- [numpy](https://numpy.org/) version 1.21.6 (BSD 3-Clause License)

# Metrics

## Tree Comparison

The comparison between phylogenetic trees (i.e., trees based on genetic data) and climatic trees involves a phylogeography step using Robinson and Foulds distance (i.e., topology distance) and Least Square distance (i.e., branch length distance).

### Least Squares Distance

$$
LS(T_1, T_2) = \sum_{i=1}^{n-1} \sum_{j=i}^{n} | \delta(i,j) - \xi(i,j) |
$$

where $T_1$ is the phylogenetic tree 1, $T_2$ is the phylogenetic tree 2, $i$ and $j$ are two species, $\delta(i,j)$ is the distance between species $i$ and species $j$ in $T_1$, $\xi(i,j)$ is the distance between species $i$ and species $j$ in $T_2$, and $n$ is the total number of species.

### Robinson-Foulds Distance

The _RF_ distance between the phylogenetic tree $T_1$ and reference tree $T_2$ is the number of non-trivial bipartitions of $$ T_1 $$ that are not in $T_2$ plus the number of non-trivial bipartitions of $T_2$ that are not in $T_1$. This distance _RF_ between $T_1$ and $T_2$ is computed by the following formula:

$$
RF(T_1,T_2) = \frac{|(Q \backslash P) \cup (P \backslash Q)|}{2n-6}
$$

where $Q$ is a set of all possible bipartitions in the phylogenetic tree (denoted $T_1$), $P$ is a set of all possible bipartitions in the reference tree (denoted $T_2$), and $n$ is the number of leaves in $T_1$ (or $T_2$). It is often relevant to normalize this distance by the maximum possible value of _RF_ (equal to $2n-6$ for two binary trees with $n$ common leaves).


# Editing Multiple Sequence Alignment Methods
Multiple Sequence Alignment (MSA) is a crucial step in bioinformatics for comparing and analyzing biological sequences. Here's an overview of some commonly used MSA methods, including pairwise alignment and popular tools like MUSCLE, CLUSTALW, and MAFFT:

- Pairwise Alignment
- MUSCLE (Multiple Sequence Comparison by Log-Expectation)
- CLUSTALW
- MAFFT (Multiple Alignment using Fast Fourier Transform)


# Conclusion

The `aPhyloGeo` pipeline represents an integrative framework, synthesizing an expansive repertoire of advanced analytical methodologies suited to diverse datasets, encompassing both genetic and climatic dimensions. By consolidating these multifaceted analyses within a unified platform, users benefit not only from simplifying the intricacies associated with navigating diverse tools but also from the assurance of heightened reproducibility in research outcomes.

In the anticipation of forthcoming advancements, `aPhyloGeo`  aims to integrate of cutting-edge functionalities. This includes the incorporation of clustering techniques predicated upon similarity derived from multiple sequence alignments, concomitant with the introduction of a more computationally efficient alignment methodology. To fortify the analytical arsenal, the integration of novel metrics, such as the Quartet metric and bipartition, aims to endow users with enhanced discernment for making nuanced decisions concerning their datasets through a comprehensive and robust assessment of genetic diversity.

By adhering resolutely adhering to the pinnacle of standards in software development, this research not only aspires to furnish an immediate solution but also endeavors to position `aPhyloGeo` as a dynamic and evolving platform. Aspiring to be a vanguard in the field of phylogeographic analysis, the pipeline is dedicated to furnishing users with a sophisticated suite of tools that seamlessly adapt to the evolving landscape of genetic research. Through these advancements, the pipeline seeks to make a substantial and enduring contribution to the corpus of scientific knowledge, elevating the benchmarks of reproducibility and usability within the broader scientific community.

# Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada, Fonds de recherche du Québec - Nature et technologie, the University of Sherbrooke grant, and the Centre de recherche en écologie de l'UdeS (CREUS). The author would like to thank the Department of Computer Science, University of Sherbrooke, Quebec, Canada for providing the necessary resources to conduct this research. The computations were performed on resources provided by Compute Canada and Compute Quebec - the National and Provincial Infrastructure for High-Performance Computing and Data Storage. The author would like to thank the students of the University of Sherbrooke and the Université du Québec à Montréal for their great contribution to the development of the software.

# References
