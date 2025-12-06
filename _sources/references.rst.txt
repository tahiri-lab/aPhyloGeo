.. _references:

References
==========

State of the Field - Advancements in Genomic Analysis
-----------------------------------------------------

The field of genomic analysis has progressed significantly in recent years, notably in the creation of tools and algorithms to explore the intricate relationship between genetic variation and environmental factors. Our algorithm for identifying sub-sequences within genes [#nadia_tahiri-proc-scipy-2022], and its subsequent application to SARS-CoV-2 data in 2023 [#nadia_tahiri-proc-scipy-2023], enhance comprehension of the genetic underpinnings of adaptation across various species and environments.

In the broader field of phylogeography, substantial methodological advancements have also occurred. Several Python packages provide functionalities pertinent to phylogeographic analysis, but often in a fragmented way. `Biopython <https://pypi.org/project/biopython>`_, a cornerstone in bioinformatics, excels at handling genetic sequences and basic phylogenetic tasks, yet falls short in integrating environmental data. `DendroPy <https://pypi.org/project/DendroPy>`_, a robust library for phylogenetic trees, aids in visualizing phylogeographic patterns but requires additional tools for comprehensive analysis. While `SciPy <https://pypi.org/project/scipy/>`_'s statistical utilities could be harnessed for custom analyses, its complexity demands a strong background in statistical programming. `GeoPandas <https://pypi.org/project/geopandas/>`_, adept at handling geospatial data, is useful for mapping genetic or environmental distributions, but lacks seamless integration with genetic data analysis tools. In summary, while powerful individual tools exist, a comprehensive and user-friendly Python package specifically designed for phylogeographic analysis remains a gap to be filled.

Statistical approaches, including generalized linear models (GLMs) and mixed models, are increasingly used to investigate the relationship between genetic variation and environmental variables. These methods enable researchers to quantify the relative influence of various factors, such as climate, geography, and demography, on observed patterns of genetic diversity.

The continuous refinement of these tools and methodologies, coupled with the growing availability of high-throughput sequencing technologies and environmental data, has opened up exciting new research avenues in evolutionary biology, ecology, and conservation. *aPhyloGeo* builds upon these advancements, providing a unified platform for integrating genetic and climatic data to address a wide array of phylogeographic questions. By bridging the gap between genomics and environmental science, *aPhyloGeo* aims to contribute to a more comprehensive understanding of the forces shaping biodiversity in a changing world.


1. **Calculation of distance between phylogenetic tree: `Least Square metric`**

   - `Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis. Models and estimation procedures. American journal of human genetics, 19(3 Pt 1), 233. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/>`_

   - `Felsenstein, J. (1997). An alternating least squares approach to inferring phylogenies from pairwise distances. Systematic biology, 46(1), 101-111. <https://pubmed.ncbi.nlm.nih.gov/11975348/>`_

   - `Makarenkov, V., & Lapointe, F. J. (2004). A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. Bioinformatics, 20(13), 2113-2121. <https://pubmed.ncbi.nlm.nih.gov/15059836/>`_

2. **Calculation of distance between phylogenetic tree: `Robinson-Foulds metric`**

   - `Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147. <https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub>`_

3. **Dataset full description: `Analysis of genetic and climatic data of SARS-CoV-2`**

   - `Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2. <https://conference.scipy.org/proceedings/scipy2022/nadia_tahiri.html>`_

   - `Li, W., & Tahiri, N. (2023). aPhyloGeo-Covid: A Web Interface for Reproducible Phylogeographic Analysis of SARS-CoV-2 Variation using Neo4j and Snakemake. <https://conference.scipy.org/proceedings/scipy2023/nadia_tahiri.html>`_

4. **Muscle5**:

   - `Edgar, RC (2021). MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping, bioRxiv 2021.06.20.449169. <https://doi.org/10.1101/2021.06.20.449169>`_

5. **Fastree**:

   - `Price, MN, Dehal, PS, and Arkin, AP (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments, PLoS ONE 5(3): e9490. <https://doi.org/10.1371/journal.pone.0009490>`_

6. **ClustalW**:

   - `Larkin, MA, Blackshields, G, Brown, NP, Chenna, R, McGettigan, PA, McWilliam, H, Valentin, F, Wallace, IM, Wilm, A, Lopez, R, Thompson, JD, Gibson, TJ, and Higgins, DG (2007). Clustal W and Clustal X version 2.0, Bioinformatics 23(21): 2947-2948. <https://doi.org/10.1093/bioinformatics/btm404>`_

7. **Mafft**:

   - `Katoh, K, Rozewicki, J, and Yamada, KD (2019). MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization, Briefings in Bioinformatics 20(4): 1160-1166. <https://doi.org/10.1093/bib/bbx108>`_

