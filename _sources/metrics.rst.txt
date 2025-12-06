.. _metrics:

Metrics
=======

Comparison Between Phylogenetic and Climatic Trees
--------------------------------------------------

The comparison between phylogenetic trees (i.e., trees based on genetic data) and climatic trees involves a phylogeography step using Robinson and Foulds distance (i.e., topology distance) and Least Square distance (i.e., branch length distance).

Least Squares Distance
-----------------------

.. math::

    LS(T_1, T_2) = \sum_{i=1}^{n-1} \sum_{j=i}^{n} \lvert \delta(i,j) - \xi(i,j) \rvert

Where :math:`T_1` is phylogenetic tree 1, :math:`T_2` is phylogenetic tree 2, :math:`i` and :math:`j` are two species, :math:`\delta(i,j)` is the distance between species :math:`i` and species :math:`j` in :math:`T_1`, :math:`\xi(i,j)` is the distance between species :math:`i` and species :math:`j` in :math:`T_2`, and :math:`n` is the total number of species.

Robinson-Foulds Distance
------------------------

The *RF* distance between the phylogenetic tree (:math:`T_1`) and reference tree (:math:`T_2`) is the number of non-trivial bipartitions of :math:`T_1` that are not in :math:`T_2` plus the number of non-trivial bipartitions of :math:`T_2` that are not in :math:`T_1`. This distance *RF* between :math:`T_1` and :math:`T_2` is computed by the following formula:

.. math::

    RF(T_1, T_2) = \frac{|(Q \backslash P) \cup (P \backslash Q)|}{2n-6},

where :math:`Q` is the set of all possible bipartitions in phylogenetic tree :math:`T_1`, :math:`P` is the set of all possible bipartitions in reference tree :math:`T_2`, and :math:`n` is the number of leaves in :math:`T_1` (or :math:`T_2`). It is often relevant to normalize this distance by the maximum possible value of *RF* (equal to :math:`2n-6` for two binary trees with :math:`n` common leaves).

