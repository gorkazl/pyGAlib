# -*- coding: utf-8 -*-
# Copyright (c) 2013, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
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

...moduleauthor:: Gorka Zamora-López <gorka@zamora-lopez.xyz>

"""
# Standard library imports
# Third party imports
import numpy as np
from numba import jit
# Local imports
from .metrics import is_directed


############################################################################
"""CONNECTIVITY AND DEGREE STATISTICS"""
def MatchingIndex(adjmatrix, directed=None, normed=False):
    """Computes the number of common neighbours for every pair of nodes.

    The matching index of two nodes i and j is the number of common
    neighbours they are connected to.

    Parameters
    ----------
    adjmatrix : ndarray of shape (N,N)
        The adjacency matrix of the network. If weighted, link weights are ignored.
    directed: None or Boolean. Optional, default: None
        If ``None``, automatically checks whether `adjmatrix` is a directed or
        an undirected graph, and continues accordingly. If ``True`` computes the
        results for both the the input and output neighbours shared by two nodes.
        If ``False``, computes the results for undirected graphs.
    normed : Boolean, optional. Default: False
        Determines whether the number or the fraction of common neighbours is
        returned. If ``False``, it computes the number of common neighbours of
        two nodes i and j. If ``True``, this is divided by the total number of
        distinct neighbours that i and j share. Thus, values range from 0 (i and
        j have no common neighbours) to 1 (i and j are connected to the same
        nodes).
        Note that explicit links i-->j and j-->i are excluded. For example:
            1 -- [2,3,4,5]
            2 -- [1,4,5,6]
            MI(1,2) = len({4,5}) / len({3,4,5,6}) = 0.5

    Returns
    -------
    mindices : ndarray of shape (N,N), or tuple with two ndarrays of shape ((N,N), (N,N)).
        The matching indices for all pairs of nodes. If ``directed = False``
        returns a single matrix and, if ``directed = True`` a tuple is returned
        containing two matrices. ``mindices[0]`` is the matrix of input pair-wise
        matching indices (the common neighbours nodes i and j receive connections
        from) and ``mindices[1]`` is the matrix for the output neighbours (number
        of common neighbours i and j project to).
        If ``directed = None``, automatically checks whether `adjmatrix` is
        directed or undirected, and return the corresponding result.

    Usage
    -----
    Except for the default ``directed=None`` case, the function does not check
    whether the input `adjmatrix` is directed or not. Therefore, for the cases
    ``directed=True`` and ``directed=False`` the function behaves as described
    above, returning either one or two matrices.
    If `adjmatrix` is a directed graph but user enters ``directed = False``,
    then the function returns a single matrix with the matching indices of the
    output neighbours (inputing the transpose would return the result for the
    input neighbours). When `adjmatrix` is an undirected graph but user enters
    ``directed = True``, then two identical matrices are returned because in a
    graph the input and output neighbours are the same.
    """
    @jit
    def _MatchingIndex_Base(boolmatrix, _normed):
        """Performs the calculation directly on the boolean matrix.
        """
        # Initialise the pair-wise matching index matrix
        N = len(boolmatrix)
        mimatrix = np.zeros((N,N), np.float64)

        for i in range(N):
            ineighbours = boolmatrix[i]
            if _normed: mimatrix[i,i] = 1.0
            else: mimatrix[i,i] = ineighbours.sum()

            for j in range(i+1,N):
                jneighbours = boolmatrix[j]
                # Intersection and union of the neighbourhoods
                mi = (ineighbours * jneighbours).sum()
                if mi == 0: continue
                mi = np.float64(mi)
                if _normed:
                    union = ineighbours + jneighbours
                    norm = union.sum()
                    # Avoid counting the explicit links i-->j and j-->i
                    if union[i] == True: norm -= 1
                    if union[j] == True: norm -= 1
                    # Normalize and save the value avoiding ZeroDivision errors
                    mi /= norm

                mimatrix[i,j] = mi
                mimatrix[j,i] = mi

        return mimatrix

    #______________________________________________________________________
    # SECURITY CHECKS AND PREPARATION
    if directed not in (None, True, False):
        raise ValueError( f"Optional parameter `directed` shall be one of (``None``, ``True``, ``False``), ``{directed}`` entered." )

    if normed not in (True, False):
        raise ValueError( f"Optional parameter `normed` must be ``True`` or ``False``, got `{normed}`" )

    # Convert to boolean and make sure self-loops are ignored
    adjmatrix = adjmatrix.astype(bool)
    np.fill_diagonal(adjmatrix, False)

    # COMPUTE THE PAIR-WISE MATRIX OF MATCHING INDICES
    if directed == None:
        directed = is_directed(adjmatrix)

    if directed == True:
        mindices_outputs = _MatchingIndex_Base(adjmatrix, normed)
        adjmatrixT = np.transpose(adjmatrix).copy()
        mindices_inputs = _MatchingIndex_Base(adjmatrixT, normed)
        mindices = (mindices_inputs, mindices_outputs)
    else:
        mindices = _MatchingIndex_Base(adjmatrix, normed)

    return mindices



###############################################################################
"""PATHS, CYCLES AND DISTANCE FUNCTIONS"""
def FloydWarshall(adjmatrix, weighted_dist=False):
    """Computes the pathlength between all pairs of nodes in a network.

    This version returns the same output as `FloydWarshall()` function in the
    metrics.py module, but runs much faster (specially for networks of N > 100)
    due to the Numba precompilation of the main loop.

    Parameters
    ----------
    adjmatrix : ndarray of shape (N,N)
        The adjacency matrix of the network.
    weighted_dist : boolean, optional
        If False, computes the binary graph distance between the nodes in the
        graph. If True, it will compute the weighted graph distance but ...

    ... WARNING !! Calculation of weighted graph distance ONLY MAKES SENSE
    WHEN the weights in `adjmatrix` represent a spatial distance between
    the nodes. For example, if `adjmatrix` is the graph of a transportation
    network and its weights encode the kilometers between the cities or
    stations (nodes) in the network.

    Returns
    -------
    distmatrix : ndarray of shape (N,N)
        The pairwise distance matrix dij of the shortest path between
        nodes i and j.
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
        distmatrix = np.where(adjmatrix == 0, np.inf, adjmatrix).astype(np.float64)
    else:
        distmatrix = np.where(adjmatrix == 0, np.inf, 1).astype(np.float64)

    # 2) RUN THE FLOYD-WARSHALL ALGORITHM USING FASTER FUNCTIONS (NUMBA)
    if is_directed(adjmatrix) == True:
        FW_Directed(distmatrix)
    else:
        FW_Undirected(distmatrix)

    return distmatrix


############################################################################
"""COMPONENTS, COMMUNITIES, K-CORES..."""


######################################################################
"""ROLES OF NODES IN NETWORKS WITH COMMUNITY (ASSORTATIVE) ORGANIZATION"""


#
