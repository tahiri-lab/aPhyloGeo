.. _sequence-similarity-algorithms:

Sequence Similarity Algorithms
==============================

Filtering Alignment Windows
---------------------------

To enhance the algorithm's efficiency, an additional step has been introduced. This step filters alignment windows, retaining only those with a dissimilarity rate surpassing the user-defined threshold specified in the parameter YAML file. These retained windows represent the most significant sequence alignments.

Edit-Based Algorithms in Bioinformatics and Phylogeny
-----------------------------------------------------

Edit-based algorithms, often referred to as distance-based algorithms, quantify the minimum number of single-character operations (insertions, deletions, or substitutions) needed to convert one sequence into another. The greater the number of operations required, the lower the similarity or greater the distance between the sequences. In the field of bioinformatics, these algorithms play a crucial role in tasks related to phylogeny and phylogeography, enabling the comparison of genetic sequences and their evolutionary relationships. They serve as a fundamental metric for various sequence similarity techniques and find widespread application in bioinformatics tasks such as DNA sequence alignment, phylogenetic analysis, and molecular evolution studies.

Hamming Distance
~~~~~~~~~~~~~~~~

The Hamming distance, denoted as :math:`H(s_1, s_2)`, measures the number of positions at which the corresponding symbols differ between two equal-length sequences. It is suitable for sequences of equal length and does not consider insertions or deletions.

.. math::

    H(s_1, s_2) = \sum_{i=1}^{n} \delta(s_{1i} \neq s_{2i})

Levenshtein Distance
~~~~~~~~~~~~~~~~~~~~

The Levenshtein distance, denoted as :math:`L(s_1, s_2)`, calculates the minimum number of single-character edits (insertions, deletions, or substitutions) needed to transform one sequence into another. It is applicable to sequences of varying lengths.


.. math::

    L(s_1, s_2) = \min\left( 
    \begin{split} 
        & L(s_1[1:m-1], s_2) + 1, \\
        & L(s_1, s_2[1:n-1]) + 1, \\
        & L(s_1[1:m-1], \\ 
        & s_2[1:n-1]) + \delta(s_{1m} \neq s_{2n}) 
    \end{split} 
    \right)


Damerau-Levenshtein Distance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Damerau-Levenshtein distance, denoted as :math:`D(s_1, s_2)`, extends the Levenshtein distance by including the transposition of adjacent characters as a valid operation. It is useful for detecting both insertions and transpositions.

.. math::

    D(s_1, s_2) = \min\left( 
    \begin{split}
        & D(s_1[1:m-1], s_2) + 1, \\
        & D(s_1, s_2[1:n-1]) + 1, \\
        & D(s_1[1:m-1], \\
        & s_2[1:n-1]) + \delta(s_{1m} \neq s_{2n}), \\
        & D(s_1[1:m-2],\\ 
        & s_2[1:n-2]) + 1  
    \end{split} 
    \right)


Jaro Similarity
~~~~~~~~~~~~~~~

Jaro similarity, denoted as :math:`J(s_1, s_2)`, measures the similarity between two strings by considering the number of matching characters and the number of transpositions. It is suitable for comparing strings with small differences and is particularly effective for short strings.

.. math::

    J(s_1, s_2) = \frac{1}{3} \left( \frac{m}{\lvert s_1 \rvert} + \frac{m}{\lvert s_2 \rvert} + \frac{m - t}{m} \right)

.. math::

    t = \frac{1}{2} \left( \text{number of transpositions} \right)

Token-Based Algorithms
----------------------

Jaccard Similarity
~~~~~~~~~~~~~~~~~~

Jaccard similarity, denoted as :math:`J(A, B)`, calculates the similarity between two sets by dividing the size of their intersection by the size of their union. In bioinformatics, it is often employed for comparing sets of elements, such as biological features.

.. math::

    J(A, B) = \frac{\lvert A \cap B \rvert}{\lvert A \cup B \rvert}

Sørensen-Dice Similarity
~~~~~~~~~~~~~~~~~~~~~~~~

Sørensen-Dice similarity, also known as Dice coefficient, is another measure for comparing the similarity between two sets. It considers the ratio of twice the size of the intersection to the sum of the sizes of the two sets.

.. math::

    S(A, B) = \frac{2 \lvert A \cap B \rvert}{\lvert A \rvert + \lvert B \rvert}

Jaro-Winkler Similarity
~~~~~~~~~~~~~~~~~~~~~~~

Jaro-Winkler similarity, denoted as :math:`JW(s_1, s_2)`, is an extension of Jaro similarity that assigns higher weights to matching prefixes. It is particularly useful when comparing strings with common prefixes.

.. math::

    JW(s_1, s_2) = J(s_1, s_2) + \ell \cdot (1 - J(s_1, s_2))

.. math::

    \ell = \text{common prefix length, maximum 4 characters}

Smith–Waterman Similarity
~~~~~~~~~~~~~~~~~~~~~~~~~

Smith–Waterman similarity, denoted as :math:`SW(s_1, s_2)`, is a local sequence alignment algorithm that identifies similar regions between sequences. It is suitable for identifying similarities in subsequences within larger sequences.

.. math::

    SW(s_1, s_2) = \max_{i,j} \left( \text{score}(s_1[i], s_2[j]) \right)

.. math::

    \text{score}(s_1[i], s_2[j]) = \max\left( 
        0, \\
        \text{score}(s_1[i-1], s_2[j-1]) + \\
        \text{similarity}(s_1[i], s_2[j]) 
    \right)

