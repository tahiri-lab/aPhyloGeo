.. _user-preferences:

Editing User Preferences
=========================

The `aPhyloGeo` software allows users to customize their preferences through a YAML configuration file. This file includes a set of parameters for easy handling. Below is an example configuration file with explanations for each parameter:

Configuration File Example
---------------------------

.. code-block:: yaml

    file_name: './datasets/example/geo.csv'
    specimen: 'id'
    dist_threshold: 60
    window_size: 200
    step_size: 100
    bootstrap_threshold: 100
    reference_gene_dir: './datasets/example'
    reference_gene_file: 'sequences.fasta'
    makeDebugFiles: True
    alignment_method: '1' # Options: 1: pairwiseAligner, 2: MUSCLE, 3: CLUSTALW, 4: MAFFT
    distance_method: '1' # Options: 1: Least-Square distance, 2: Robinson-Foulds distance, 3: Euclidean distance (DendroPY)
    fit_method: '1' # Options: 1: Wider Fit by elongating with Gap (starAlignment), 2: Narrow-fit prevent elongation with gap when possible
    tree_type: '1' # Options: 1: BioPython consensus tree, 2: FastTree application
    rate_similarity: 90
    method_similarity: '1' # Options: 1: Hamming distance, 2: Levenshtein distance, 3: Damerau-Levenshtein distance, 4: Jaro similarity, 5: Jaro-Winkler similarity, 6: Smith–Waterman similarity, 7: Jaccard similarity, 8: Sørensen-Dice similarity
    preprocessing_genetic: 1                # Enable genetic preprocessing (1 = yes, 0 = no)
    preprocessing_climatic: 1               # Enable climatic preprocessing (1 = yes, 0 = no)
    preprocessing_threshold_genetic: 0.2    # Proportion of gaps allowed per column in genetic alignments
    preprocessing_threshold_climatic: 0.7   # Variance threshold for filtering climatic features

User Preferences Options
-------------------------

- **File Name**: Path to the input data file (**`./datasets/example/geo.csv`** in the example).

  
- **Specimen**: Identifier for the specimens in the dataset (**`id`** in the example).

  
- **Bootstrap Threshold**: Number of replicates threshold to be generated for each sub-MSA.

  
- **Distance Threshold**: Distance threshold between genetic tree and climatic tree for each sub-MSA.

  
- **Window Size**: Size of the sliding window.

  
- **Step Size**: Sliding window advancement step.


  
- **Data Names**: List of newick file names for each dataset.

  
- **Reference Gene Directory**: Directory containing reference gene data (**`./datasets/example`** in the example).

  
- **Reference Gene File**: File containing reference gene sequences (**`sequences.fasta`** in the example).

  
- **Make Debug Files**: Option to generate debug files (**True** or **False**).

  
- **Alignment Method**: Algorithm selection for sequence alignment (Options: **`1: pairwiseAligner, 2: MUSCLE, 3: CLUSTALW, 4: MAFFT`** in the example).

  
- **Distance Method**: Distance selection (Options: **`1: Least-Square distance, 2: Robinson-Foulds distance, 3: Euclidean distance (DendroPY)`** in the example).

  
- **Fit Method**: Gap selection elongation (Options: **`1: Wider Fit by elongating with Gap (starAlignment), 2: Narrow-fit prevent elongation with gap when possible`** in the example).

  
- **Tree Inference Method**: The choice of inference methods (Options: **`1: BioPython consensus tree, 2: FastTree application`** in the example).

  
- **Rate Similarity**: The rate similarity between sequences to reduce and remove the sub-MSA with a high value of similarity.

  
- **Method Similarity**: The choice of similarity methods (Options: **`1: Hamming distance, 2: Levenshtein distance, 3: Damerau-Levenshtein distance, 4: Jaro similarity, 5: Jaro-Winkler similarity, 6: Smith–Waterman similarity, 7: Jaccard similarity, 8: Sørensen-Dice similarity`** in the example).

- **Genetic Preprocessing**: Enable or disable filtering of alignment columns with gaps (**`preprocessing_genetic: 1`** to enable, **`0`** to disable).

- **Climatic Preprocessing**: Enable or disable variance-based filtering of climatic features (**`preprocessing_climatic: 1`** to enable, **`0`** to disable).

- **Genetic Preprocessing Threshold**: Maximum allowed proportion of gaps per column in the alignment (**`preprocessing_threshold_genetic: 0.2`** in the example).

- **Climatic Preprocessing Threshold**: Minimum variance threshold to retain a climatic feature (**`preprocessing_threshold_climatic: 0.7`** in the example).

To use the following alignement methods, **MUSCLE**, **CLUSTALW**, and **MAFFT**, please ensure to follow the installation instructions provided in the `Alignment Dependencies Installation <alignment_dependencies.html>`_ section.