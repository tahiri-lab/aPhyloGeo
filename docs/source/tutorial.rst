.. _tutorial:

Tutorial
========

Description
-----------

Cumacea (Crustacea: Peracarida) from the deep North Atlantic to the Arctic Ocean
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The study area was located in a northern region of the North Atlantic, including the Icelandic Sea, the Denmark Strait, and the Norwegian Sea. The specimens examined were collected as part of the IceAGE project (Icelandic marine Animals: Genetic and Ecology; Cruise ship M85/3 in 2011), which studied the deep continental slopes and abyssal waters around Iceland `Meißner et al., 2018 <https://doi.org/10.1007/s12526-018-0884-7>`_. The sampling period for the included specimens was from August 30 to September 22, 2011, and they were collected at depths ranging from 316 to 2568 m. Information on the sampling plan, sample processing, DNA extraction steps, PCR amplification, sequencing, and extracted and aligned DNA sequences is available in the `Uhlir et al., 2021 <https://doi.org/10.7717/peerj.12379>`_ article. Refer to the `example Cumacea FASTA file <https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/Cumacea/Cumacea.fasta>`_ and `example Cumacea CSV file <https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/Cumacea/Cumacea.csv>`_ for guidance.

Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a previous study (Koshkarov et al., 2022), a total of 38 distinct genetic lineages were identified across the species' range. For this specific analysis, we focused on 5 of these lineages, selected for their pronounced regional characteristics and relevance to our research. Based on location information, complete nucleotide sequencing data for these 5 lineages was collected from the `NCBI Virus website <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/>`_. In the case of the availability of multiple sequencing results for the same lineage in the same country, we selected the sequence whose collection date was closest to the earliest date presented. If there are several sequencing results for the same country on the same date, the sequence with the least number of ambiguous characters (N per nucleotide) is selected.

Although the selection of samples was based on the phylogenetic cluster of lineage and transmission, most of the sites involved represent different meteorological conditions. The 5 samples involved temperatures ranging from -4°C to 32.6°C, with an average temperature of 15.3°C. The specific humidity ranged from 2.9 g/kg to 19.2 g/kg with an average of 8.3 g/kg. The variability of wind speed and all-sky surface shortwave downward irradiance was relatively small across samples compared to other parameters. The wind speed ranged from 0.7 m/s to 9.3 m/s with an average of 4.0 m/s, and all-sky surface shortwave downward irradiance ranged from 0.8 kW-hr/m²/day to 8.6 kW-hr/m²/day with an average of 4.5 kW-hr/m²/day. In contrast to the other parameters, 75% of the cities involved receive less than 2.2 mm of precipitation per day, and only 5 cities have more than 5 mm of precipitation per day. The minimum precipitation is 0 mm/day, the maximum precipitation is 12 mm/day, and the average value is 2.1 mm/day. Refer to the `example SARS-CoV-2 FASTA file <https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/example/sequences.fasta>`_ and `example SARS-CoV-2 CSV file <https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/example/geo.csv>`_ for guidance.

Input
------

The algorithm takes two files as input with the following definitions:

- 🧬 **Genetic file** with FASTA extension. The first file or set of files will contain the genetic sequence information of the species sets selected for the study. The name of the file must allow you to know the name of the gene. It is therefore strongly recommended to follow the following nomenclature: ``gene_name.fasta``. It should contain genetic variants (e.g., SNPs) and their associated metadata (e.g., sample IDs, location information).
- ⛅ **Climatic file** with CSV extension (Comma-Separated Values). The second file will contain the habitat information for the species sets selected for the study. Each row will represent the species identifier and each column will represent a climate condition. It should include relevant climatic variables (e.g., temperature, precipitation) for each geographic location represented in your genetic data and must be clearly labeled to match the expected format.

Preparing Your Data
~~~~~~~~~~~~~~~~~~~

Include relevant climatic variables (e.g., temperature, precipitation) for each geographic location represented in your genetic data. Column headers must be clearly labeled to match the expected format. Refer to the example files in the datasets directory for guidance.

Output
------

The algorithm will return a CSV file that contains information from all relevant MSAs (see the Workflow section for more details). The sliding windows of interest are those with interesting bootstrap support (i.e., indicating the robustness of the tree) and high similarity to the climate condition in question (i.e., based on the ``RF``, ``RFnorm``, ``LS``, and ``Euclidean`` values). They will indicate, among other things, the name of the gene, the position of the beginning and end of the sliding window, the average bootstrap value, the LS value, and finally the climatic condition for which this genetic zone would explain the adaptation of the species in a given environment.

To sum up, aPhyloGeo generates an ``output.csv`` file containing analysis results. Additional visualizations (e.g., maps, plots) may be generated based on your configuration.

Prerequisites
-------------

System Requirements
~~~~~~~~~~~~~~~~~~~

* **Operating System:** Windows, macOS, or Linux.
* **Python:** `Python <https://www.python.org>`_ 3.8 or higher.

Key Features of aPhyloGeo
~~~~~~~~~~~~~~~~~~~~~~~~~

* **Multi-Platform:** Works seamlessly on Windows, macOS, and Linux.
* **Flexible Analysis:** Supports various phylogeographic analyses, including:
    - Identifying genetic lineages and their geographic origins
    - Assessing the impact of climate on genetic diversity
    - Visualizing genetic and geographic relationships
* **Customizable:** Tailor analyses using a configuration file to fit your specific research questions.
* **Open Source:** Freely available and encourages contributions from the research community.

Before you begin, ensure you have the following installed:

.. code-block:: bash

    pip install pandas aphylogeo

