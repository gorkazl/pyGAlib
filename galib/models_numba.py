# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
FASTER SYNTHETIC NETWORK GENERATORS
===================================

This module duplicates some functions in the module 'gamodels.py'
and significantly accelerates their performance. This is achieved thanks to
the Numba package (http://numba.pydata.org) which compiles annotated Python
and NumPy code to LLVM on run-time. Therefore, the use of this module
requires Numba to be installed.

Due to the run-time approach of Numba, the original functions in 'galib.py'
will run faster in very small networks and Numba-based functions take the
lead as network size increases. The precise size limit varies from function
to function and I will report approximate valuesfor each of them. In general,
for any network of N > 100 nodes Numba-based functions run faster.

RANDOM NETWORK GENERATORS
-------------------------
RandomGraph_Numba
    Generates random graphs with N nodes and L links.


...moduleauthor:: Gorka Zamora-Lopez <galib@zamora-lopez.xyz>

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import numpy.random
from numba import jit


############################################################################
"""RANDOM NETWORK GENERATORS"""
@jit
def RandomGraph_Numba(N, L, directed=False, selfloops=False):
    """Generates random graphs with N nodes and L links.

    Similar to an Erdos-Renyi (ER) graph with probability p = rho, where
    rho is the density of links. In ER graphs the total number of links
    varies in different realizations what is unsuitable to compare empirical
    networks with their random counterparts. RandomGraph() allows to create
    random graphs that always have the same number of links. The precise
    formula for rho depends on the options given. For an undirected graph
    with no self-loops allowed, rho = 1/2 * L / (N*(N-1)).

    Parameters
    ----------
    N : integer
        The size of the network (number of nodes).
    L : integer
        Number of links of the resulting random network.
    directed : Boolean
        True if a directed graph is desired. False, for an undirected graph.
    selfloops: Boolean
        True if self-loops are allowed, False otherwise.

    Returns
    -------
    adjmatrix : ndarray of rank-2, size NxN and dtype = int.
        The adjacency matrix of the generated random graph.

    Notes
    -----
    Make sure to specify the right number of links for the type of graph
    desired. Keep in mind that directed graphs have twice as many links
    (arcs) as graphs (edges).

    See Also
    --------
    ErdosRenyiGraph : Random graphs with given link probability.
    """
    # 0) SECURITY CHECKS. Make sure L is not too large.
    if directed:
        if selfloops:
            maxL = N**2
            if L > maxL:
                raise ValueError("L out of bounds, max(L) = N**2 = %d" %maxL)
        else:
            maxL = N*(N-1)
            if L > maxL:
                raise ValueError("L out of bounds, max(L) = N*(N-1) = %d" %maxL)
    else:
        if selfloops:
            maxL = 0.5*N*(N+1)
            if L > maxL:
                raise ValueError("L out of bounds, max(L) = 1/2*N*(N+1) = %d" %maxL)
        else:
            maxL = 0.5*N*(N-1)
            if L > maxL:
                raise ValueError("L out of bounds. For the options given, max(L) = 1/2*N*(N-1) = %d" %maxL)

    # 1) INITIATE THE MATRIX AND HELPERS
    adjmatrix = np.zeros((N,N), np.uint8)
    counter = 0

    # 2) GENERATE THE MATRIX
    while counter < L:
        # 2.1) Pick up two nodes at random
        source = int(N * numpy.random.rand())
        target = int(N * numpy.random.rand())

        # 2.2) Check if they can be linked, otherwise look for another pair
        if adjmatrix[source,target] == 1: continue
        if not selfloops and source == target: continue

        # 2.3) If the nodes are linkable, place the link
        adjmatrix[source,target] = 1
        if not directed:
            adjmatrix[target,source] = 1

        counter += 1

    return adjmatrix


############################################################################
"""NETWORK REWIRING ALGORITHMS"""


############################################################################
"""MODULAR AND HIERARCHICAL NETWORK MODELS"""


##
