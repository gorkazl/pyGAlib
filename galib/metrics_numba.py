# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
"""
FASTER GRAPH ANALYSIS DESCRIPTORS
=================================

This module duplicates the slowest functions of main module 'metrics.py'
and significantly accelerates their performance. This is achieved thanks to
the Numba package (http://numba.pydata.org) which compiles annotated Python
and NumPy code to LLVM on run-time. Therefore, the use of this module
requires Numba to be installed.

Due to the run-time approach of Numba, the original functions in 'metrics.py'
will run faster in very small networks and Numba-based functions take the
lead as network size increases. The precise size limit varies from function
to function and I will report approximate valuesfor each of them. In general,
for any network of N > 100 nodes Numba-based functions run faster.

BASIC CONNECTIVITY DESCRIPTORS
------------------------------
MatchingIndex_Numba
    Computes the number of common neighbours of every pair of nodes.

PATHS AND GRAPH DISTANCE FUNCTIONS
----------------------------------
FloydWarshall_Numba
    Computes the pathlength between all pairs of nodes in a network.


...moduleauthor:: Gorka Zamora-Lopez <galib@zamora-lopez.xyz>

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numba import jit
from . import metrics

# __all__ = ['MatchingIndex_Numba', 'FloydWarshall_Numba']


############################################################################
"""CONNECTIVITY AND DEGREE STATISTICS"""
@jit
def MatchingIndex_Numba(adjmatrix, normed=True):
    """Computes the number of common neighbours of every pair of nodes.

    The matching index of two nodes i and j is the number of common
    neighbours they are linked with.

    Returns same result as MatchingIndex() in metrics.py but uses an ndarray-based
    algorithm and the Numba package to accelerate the calculations. Recommended
    only for networks of size N > 250 nodes.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.

    normed : Boolean, optional
        If 'normed=False', returns the number of common neighbours
        of two nodes i and j. If 'normed=True', the number of common
        neighbours is divided by the total number of different nodes they
        are linked with, so values range from 0 (no common neighbours)
        to 1 (all neighbours are common to i and j). Explicit links i-->j
        and j-->i are excluded because they do not contribute to the
        matching. For example:
            1 --> [2,3,4,5]
            2 --> [4,5,6]
            MI(1,2) = len([4,5])/len([3,4,5,6]) = 0.5

    Returns
    -------
    MImatrix : ndarray of rank-2 of ndtype 'float64'
        A matrix containing the matching index for all pairs of nodes.

    Notes
    -----
    - The function accepts weighted networks but it ignores the weights.
    - If adjmatrix is directed, calling MatchingIndex(adjmatrix) computes the
    matching of the output neighbours. Passing the transpose adjmatrix.T to
    the function computes the matching of the input neighbours.
    """
    N = len(adjmatrix)
    # Convert the adjacency matrix into a boolean matrix
    adjmatrix = adjmatrix.astype(np.bool)

    MImatrix = np.identity(N, np.float64)

    for i in range(N):
        for j in range(i,N):

            # Number of common neighbours (intersection of neighbourhoods)
            mi = (adjmatrix[i] * adjmatrix[j]).sum()

            if normed:
                # Size of the union of the neighbourhoods
                norm = (adjmatrix[i] + adjmatrix[j]).sum()
                # Avoid counting the explicit links i-->j and j-->i
                if adjmatrix[i,j]: norm -= 1
                if adjmatrix[j,i]: norm -= 1

                # Normalize and save the value avoiding ZeroDivision errors
                if norm > 0:
                    mi = np.float64(mi) / norm
                    MImatrix[i,j] = mi
                    MImatrix[j,i] = mi
            else:
                # Save the value
                mi = np.float64(mi)
                MImatrix[i,j] = mi
                MImatrix[j,i] = mi

    return MImatrix


###############################################################################
"""PATHS, CYCLES AND DISTANCE FUNCTIONS"""
def FloydWarshall_Numba(adjmatrix, weighted_dist=False):
    """Computes the pathlength between all pairs of nodes in a network.

    WARNING! This version returns the same output as 'FloydWarshall()'
        function in main metrics.py module but runs much faster (for networks
        of N > 100). It requires package Numba to be installed.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    weighted_dist : boolean, optional
        If True, if the graph distances are computed considering the weights
        of the links. False, otherwise. If 'adjmatrix' is a weighted
        network but'weighted_dist = False', the weights of the links are
        ignored.

    Returns
    -------
    distmatrix : ndarray of rank-2
        The pairwise distance matrix dij of the shortest pathlength between
        nodes i and j.

    See Also
    --------
    FloydWarshall : Computes the pathlength between all pairs of nodes.
    """
    # 0) DEFINE THE CORE OF THE FW ALGORITHM, ACCELARATED BY 'Numba'
    @jit
    def FW_Undirected(distmatrix):
        """The Floyd-Warshall algorithm for undirected networks
        """
        N = len(distmatrix)
        for k in range(N):
            for i in range(N):
                for j in range(i,N):
                    d = distmatrix[i,k] + distmatrix[k,j]
                    if distmatrix[i,j] > d:
                        distmatrix[i,j] = d
                        distmatrix[j,i] = d

    @jit
    def FW_Directed(distmatrix):
        """The Floyd-Warshall algorithm for directed networks
        """
        N = len(distmatrix)
        for k in range(N):
            for i in range(N):
                for j in range(N):
                    d = distmatrix[i,k] + distmatrix[k,j]
                    if distmatrix[i,j] > d:
                        distmatrix[i,j] = d

    ########################################################################
    # 1) PREPARE FOR THE CALCULATIONS
    # 1.1) Initialize the distance matrix
    if weighted_dist:
        distmatrix = np.where(adjmatrix == 0, np.inf, adjmatrix)
    else:
        distmatrix = np.where(adjmatrix == 0, np.inf, 1)

    # 1.2) Find out whether the network is directed or undirected
    recip = metrics.Reciprocity(adjmatrix)

    # 2) RUN THE FLOYD-WARSHALL ALGORITHM USING FASTER FUNCTIONS (NUMBA)
    if recip==1.0:
        FW_Undirected(distmatrix)
    else:
        FW_Directed(distmatrix)

    return distmatrix


############################################################################
"""COMPONENTS, COMMUNITIES, K-CORES..."""


######################################################################
"""ROLES OF NODES IN NETWORKS WITH COMMUNITY (ASSORTATIVE) ORGANIZATION"""


#
