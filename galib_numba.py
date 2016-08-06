"""
=================================
FASTER GRAPH ANALYSIS DESCRIPTORS
=================================

This module duplicates the slowest functions of main module 'galib.py'
and significantly accelerates their performance. This is achieved thanks to
the Numba package (http://numba.pydata.org) which compiles annotated Python
and NumPy code to LLVM on run-time. Therefore, the use of this module
requires Numba to be installed.

Due to the run-time approach of Numba, the original functions in 'galib.py'
will run faster in very small networks and Numba-based functions take the
lead as network size increases. The precise size limit varies from function
to function and I will report approximate valuesfor each of them. In general,
for any network of N > 100 nodes Numba-based functions run faster.

BASIC CONNECTIVITY DESCRIPTORS
==============================
None yet.

PATHS AND GRAPH DISTANCE FUNCTIONS
==================================
FloydWarshall
    Computes the pathlength between all pairs of nodes in a network.

COMMUNITIES, COMPONENTS, K-CORES, ...
=====================================
None yet.

ROLES OF NODES IN NETWORKS WITH MODULAR ORGANIZATION
====================================================
None yet.
"""

__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@Zamora-Lopez.xyz"
__copyright__ = "Copyright 2013-2016"
__license__ = "GPL"
__update__="05/08/2016"

import numpy as np
from numba import autojit
import gatools
import galib


############################################################################
"""CONNECTIVITY AND DEGREE STATISTICS"""

###############################################################################
"""PATHS, CYCLES AND DISTANCE FUNCTIONS"""
def FloydWarshall_Numba(adjmatrix, weighted_dist=False):
    """Computes the pathlength between all pairs of nodes in a network.

    WARNING! This version returns the same output as 'FloydWarshall()'
        function in main galib.py module but runs much faster (for networks
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
    @autojit
    def FW_Undirected(distmatrix):
        """The Floyd-Warshall algorithm for undirected networks
        """
        N = len(distmatrix)
        for k in xrange(N):
            for i in xrange(N):
                for j in xrange(i,N):
                    d = distmatrix[i,k] + distmatrix[k,j]
                    if distmatrix[i,j] > d:
                        distmatrix[i,j] = d
                        distmatrix[j,i] = d

    @autojit
    def FW_Directed(distmatrix):
        """The Floyd-Warshall algorithm for directed networks
        """
        N = len(distmatrix)
        for k in xrange(N):
            for i in xrange(N):
                for j in xrange(N):
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
    recip = galib.Reciprocity(adjmatrix)

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
