# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
GRAPH ANALYSIS DESCRIPTORS
==========================

This module contains functions to compute many graph descriptors of a
network. Many functions compute also the results both for graphs and
digraphs. Most of functions accept weighted networks as input but they
ignore the weights unless explicitely specified. Support for weighted
measures will be added to GAlib in future releases.

BASIC CONNECTIVITY DESCRIPTORS
------------------------------
Density
    Returns the density of links in a network.
Degree
    Computes the number of neighbours of every node.
Intensity
    The total strength of a node in a weighted network.
Reciprocity
    Computes the fraction of reciprocal links to total number of links.
ReciprocalDegree
    Returns the reciprocal degree and excess degrees of every nodes.
AvNeighboursDegree
    Average neighbours' degree of nodes with given degree k, for all k.
Clustering
    Returns the clustering coefficient and the local clustering of every node.
k_Density
    Computes the density of subnetworks with degree >= k', for all k' = 0 to kmax.
RichClub
    Identifies the subset of hubs with dense interconnectivity.
MatchingIndex
    Computes the number of common neighbours of every pair of nodes.

PATHS AND GRAPH DISTANCE FUNCTIONS
----------------------------------
FloydWarshall
    Computes the pathlength between all pairs of nodes in a network.
PathsAllinOne
    Returns pathlength and betweenness. Finds all shortest paths and cycles.
AllShortestPaths
    Finds all the shortest paths between two nodes.

COMMUNITIES, COMPONENTS, K-CORES, ...
-------------------------------------
AssortativityMatrix
    Returns the assortativity matrix of network given a partition of nodes.
ConnectedComponents
    Finds all the connected components in a network out of a distance matrix.
Modularity
    Computes the Newman modularity given a partition of nodes.
K_Core
    Finds the K-core of a network with degree k >= kmin.
K_Shells
    Returns the K-shells of a network for all k from kmin to kmax.

ROLES OF NODES IN NETWORKS WITH MODULAR ORGANIZATION
----------------------------------------------------
GlobalHubness
    Computes the global hubness of all nodes in a network.
LocalHubness
    Given a partition, computes the local hubness of all nodes.
ParticipationMatrix
    Given a partition of the network, it returns the participation matrix.
ParticipationVectors
    Computes the probability of nodes to belong to every community.
NodeParticipation
    Participation index of every node given a partition of the network.
NodeDispersion
    Dispersion index of every node given a partition of the network.
NodeRoles
    Computes all four parameters to characterise the roles of nodes.
ParticipationIndex_GA
   Returns the participation index as defined by Guimera & Amaral.
Hubness_GA
    Returns the within-module degree defined by Guimera & Amaral.


...moduleauthor:: Gorka Zamora-Lopez <galib@zamora-lopez.xyz>

"""
from __future__ import division, print_function, absolute_import

import types
import numpy as np

from . import tools


############################################################################
"""CONNECTIVITY AND DEGREE STATISTICS"""
def Density(adjmatrix):
    """Returns the density of links in a network.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.

    Returns
    -------
    A scalar value between 0 and 1.
    """
    N = len(adjmatrix)
    L = adjmatrix.astype('bool').sum()

    if adjmatrix.trace():
        return float(L) / N**2
    else:
        return float(L) / (N*(N-1))

def Degree(adjmatrix, directed=False):
    """Computes the number of neighbours of every node.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    directed: Boolean, optional
        True if the network is directed, False otherwise.

    Returns
    -------
    If 'directed=False':
    degarray : ndarray
        An array with the degree of every node in the undirected network.
    If 'directed=True'
    degarray : A tuple containing two ndarrays, indegarray and outdegarray.
        indegarray is the input degree of the every node and outdegarray is
        is the output degree of every node.

    See Also
    --------
    Intensity : Computes the weighted degree of networks.
    ReciprocalDegree : Reciprocity of every node and excess degrees.
    """
    N = len(adjmatrix)
    adjmatrix = adjmatrix.astype('bool')

    if directed:
        indegree = adjmatrix.sum(axis=0)
        outdegree = adjmatrix.sum(axis=1)
        return indegree, outdegree

    else:
        degree = adjmatrix.sum(axis=1)
        return degree

def Intensity(adjmatrix, directed=False):
    """Computes the total strength of a node in a weighted network.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    directed: Boolean. optional
        True if the network is directed, False otherwise.

    Returns
    -------
    If 'directed=False'
    intensity : ndarray.
        The weighted degree of every node in the undirected network.
    If 'directed=True'
    intensity : tuple containing two ndarrays, inintensity and outintensity
        inintensity is the weighted input degree of every node in the
        directed network, outintensity is the output degree of every node
        in the directed network.

    See Also
    --------
    Degree : Computes the degree of every node in the network.
    """
    N = len(adjmatrix)

    if directed:
        inintensity = adjmatrix.sum(axis=0)
        outintensity = adjmatrix.sum(axis=1)
        return inintensity, outintensity

    else:
        intensity = adjmatrix.sum(axis=0)
        return intensity

def Reciprocity(adjmatrix):
    """Computes the fraction of reciprocal links to total number of links.

    Both weighted and unweighted input matrices are permitted. Weights
    are ignored for the calculation.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.

    Returns
    -------
    reciprocity : float
        A scalar value between 0 (for acyclic directed networks) and 1 (for
        fully reciprocal).
    """
    # 0) PREPARE FOR COMPUTATIONS
    adjmatrix = adjmatrix.astype('bool')

    # 1) COMPUTE THE RECIPROCITY
    L = adjmatrix.sum()
    if L == 0:
        reciprocity = 0
    else:
        # Find the assymmetric links
        # Rest = np.abs(adjmatrix - adjmatrix.T)
        Rest = np.abs(adjmatrix ^ adjmatrix.T)
        Lsingle = 0.5*Rest.sum()
        reciprocity = np.float(L-Lsingle) / L

    return reciprocity

def ReciprocalDegree(adjmatrix, normed=False):
    """Returns the reciprocal degree and excess degrees of every nodes.

    The reciprocal degree, kr, of a node i is the number of neighbours with
    which i makes reciprocal connections. The excess input degree, k-,
    is the number of links j-->i for which no reciprocal j<--i exist.
    k- = in-k - k_r. The excess output degree is the number of links i-->j for
    which no reciprocal i<--j link exists. k+ = out-k - kr.
    In case of undireceted networks, kr = degree, k+ = k- = 0.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    normed: Boolean. optional
        True if normalized output is desired, False otherwised.

    Returns
    --------
    A tuple containing three arrays (kr, k-, k+). The arrays give the kr, k-
    and k+ values for all the N nodes.
    If the option 'normed=True' then kr, k- and k+ are given in fractions.
        nkr = 2 * kr / (in-k + out-k)
        nk- = k- / in-k
        nk+ = k+ / out-k

    See Also
    --------
    Reciprocity : Computes the reciprocity of a directed network.

    Notes
    -----
    Both weighted and unweighted adjacency matrices are permitted, but the
    computation ignores the weights of the arcs in the case of weighted
    adjacency matrices.
    """
    # Check whether the matrix is binary of weighted
    if adjmatrix.max() != 1:
        adjmatrix = np.where(adjmatrix != 0, 1, 0)

    # Compute the input and output degrees
    indegree, outdegree = Degree(adjmatrix, True)

    # Find the symmetric and assymmetric links
    rest = abs(adjmatrix - adjmatrix.T)
    recipadjmatrix = (adjmatrix + adjmatrix.T - rest)/2
    del rest

    if normed:
        # Normalize the reciprocal degrees of the nodes
        recipdegree = np.add.reduce(recipadjmatrix)
        degminus = indegree - recipdegree
        degplus = outdegree - recipdegree

        return 2.0 * recipdegree.astype(np.float) / (indegree+outdegree), \
               degminus.astype(np.float)/indegree, \
               degplus.astype(np.float)/outdegree

    else:
        # The reciprocal degree of the nodes
        recipdegree = np.add.reduce(recipadjmatrix)
        degminus = indegree - recipdegree
        degplus = outdegree - recipdegree

        return recipdegree, degminus, degplus

def AvNeighboursDegree(adjmatrix, knntype='undirected', fulloutput=False):
    """Average neighbours' degree of nodes with given degree k, for all k.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    knntype : string
        Defines the class of degree-degree correlation desired.
        -'undirected' Use this only if 'adjmatrix' is an undirected adjacency
        matrix. Otherwise an error will be raised.
        -'outout' computes the average output degree <out-k'> of the output
        neighbours of nodes with output degree out-k.
        -'outin' computes the average input degree <in-k'> of the output
        neighbours of nodes with output degree <out-k>
        -'inin' computes the average input degree <in-k'> of the input
        neighbours of nodes with input degree in-k.
        -'inout' computes the average output degree <outk-k'> of the input
        neighbours of nodes with input degree in-k.
        -'average' computes the average mean degree 1/2*(in-k' + out-k')
        of nodes with mean degree 1/2*(in-k + out-k).
    fulloutput : boolean
        If True, it also returns 'alldata'.

    Returns
    --------
    AvKnn : ndarray of dtype 'float'. AvKnn has three rows:
            Row-0: the (output) degree, k'
            Row-1: average (input) degree of all neighbours of nodes with
            degree k = k'
            Row-2: the standard deviation of the average in row-1
    alldata : ndarray with two rows of dtype float or integer.
        The list of all node's degree k and all their neighbours' degree k'.
            Row-0: the degree k of a node
            Row-1: the degree k' of one of its neighbours

    Notes
    -----
    1) In order to compute the standard deviation the algorithm first pools
    together all neighbours' degree of all nodes with degree k and then
    computes the average and std values.

    2) The function accepts weighted adjacency matrices but it ignores the
    weights of the links. The function works with both directed and with
    undirected networks.

    3) In the case of directed networks there are four classes of degre-degree
    correlations. This can be confusing. Make sure to understand which of
    the four correlation classes you are asking the function to compute.

    4) 'alldata' is provided only for plotting purposes, such that the average
    neighbours' degree can be plotted with all the datapoints in the background.
    See documentation and the script 'Example_NeighboursDegree.py'.
    """

    # 0) Security checks and prepare data for calculations
    keylist = ('undirected', 'outin',  'outout', 'inout', 'inin', 'average')
    if knntype not in keylist:
        raise KeyError("Enter a valid knntype:", keylist)

    N = len(adjmatrix)

    if knntype == 'undirected':
        if Reciprocity(adjmatrix) < 1.0:
            raise TypeError("Option 'undirected' requires an undirected adjacency matrix")
        indegree, outdegree = Degree(adjmatrix, True)
    elif knntype == 'outin':
        indegree, outdegree = Degree(adjmatrix, True)
    elif knntype == 'outout':
        indegree, outdegree = Degree(adjmatrix, True)
        indegree = outdegree
    elif knntype == 'inout':
        adjmatrix = adjmatrix.T
        indegree, outdegree = Degree(adjmatrix, True)
    elif knntype == 'inin':
        adjmatrix = adjmatrix.T
        indegree, outdegree = Degree(adjmatrix, True)
        indegree = outdegree
    elif knntype == 'average':
        indeg, outdeg = Degree(adjmatrix, True)
        avdeg = 0.5 * (indeg + outdeg)
        indegree = avdeg; outdegree = avdeg
        adjmatrix = SymmetriseMatrixAverage(adjmatrix)

    # 1) Find the degrees k of all the neighbours of nodes with given degree k'
    kdict = {}
    for i, kout in enumerate(outdegree):
        if kout not in kdict: kdict[kout] = []
        neighbours = adjmatrix[i].nonzero()[0]
        for node in neighbours:
            kinneigh = indegree[node]
            kdict[kout].append(kinneigh)

    # 2) Compute the av. degree for all neighbours of nodes with degree k'
    klist = np.sort(list(kdict.keys()))
    # klist = (kdict.keys()).sort()
    AvKnn = np.zeros((3,len(klist)), np.float)
    for count, k in enumerate(klist):
        dummy = np.array(kdict[k])
        avk = dummy.mean()
        devk = dummy.std()

        AvKnn[0,count] = k
        AvKnn[1,count] = avk
        AvKnn[2,count] = devk

    if fulloutput:
        L = indegree.sum()
        # Convert the kdict into a 2D array for plotting purposes
        if knntype == 'Mixed': neighkarray = np.zeros((2,L),np.float)
        else: neighkarray = np.zeros((2,L),np.int)

        counter = 0
        for k in kdict:
            for node in kdict[k]:
                neighkarray[0,counter] = k
                neighkarray[1,counter] = indegree[node]
                counter += 1
        return AvKnn, neighkarray

    else:
        return AvKnn

def Clustering(adjmatrix, checkdirected=True):
    """Returns the clustering coefficient and the local clustering of every node.

    Directed networks are NOT accepted by the function. Weighted networks are
    accepted but the weights of the links will be ignored.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    checkdirected : boolean (optional)
        If 'True', it calls Reciprocity() to make sure adjmatrix is undirected.
        If 'False', it skips the check and runs a bit faster. Only recommended
        for large networks and in exceptional cases when you exactly know the
        input adjacency matrix is undirected.

    Returns
    --------
    A tuple containing two elements:
    C : Scalar between 0 and 1.
        The clustering coefficient of the network (2 * Ntriangles / Ntriads).
    ClustNodes : ndarray of dtype=float.
        An array of length N with the clustering of every node.

    Notes
    -----
    The clustering coefficient C is not the same as the average of
    the clusterings of the individual nodes <Ci>!! Here the array ClustNodes
    is provided because it costs no additional resources and the data is
    useful for an further statistical analysis.
    """
    N = len(adjmatrix)
    adjmatrix = np.where(adjmatrix,1,0).astype(np.float32)

    # 0) SECURITY CHECKS
    if checkdirected:
        if Reciprocity(adjmatrix) < 1.0:
            raise TypeError("Please introduce an undirected adjacency matrix.")

    # Remove diagonal entries, in case there is any self-loop
    adjmatrix[np.diag_indices(len(adjmatrix))] = 0

    # 1) COMPUTE THE NUMBER OF TRIANGLES EACH NODE PARTICIPATES IN
    pow2matrix = np.dot(adjmatrix,adjmatrix)
    ntriangles = np.zeros(N, np.float64)
    for i in range(N):
        ntriangles[i] = (adjmatrix[i] * pow2matrix[i]).sum()

    # 2) COMPUTE THE NUMBER OF DIADS EACH NODE PARTICIPATES IN
    deg = adjmatrix.sum(axis=1)
    ndiads = deg*(deg-1)

    # 3) COMPUTE THE COEFFICIENT AND THE CLUSTERING OF EACH NODE
    # The usual multiplication x3 is not needed because by computing
    # the matrix power A**3, we have already included it.
    # Avoid Error for very sparse networks with no diads
    Ndiads = ndiads.sum()
    if Ndiads == 0:
        coefficient = 0
        cnodes = np.zeros(N,float)
    # Compute clustering if at least one triad has been found
    else:
        coefficient = ntriangles.sum() / Ndiads
        if 0 in ndiads:
            ndiads = np.where(ndiads==0,1,ndiads)
        cnodes = ntriangles / ndiads

    return coefficient, cnodes

def k_Density(adjmatrix, rctype='undirected'):
    """Computes the density of subnetworks with degree >= k', for all k' = 0 to kmax.

    The k-density (phi(k)) is the density of the subnetwork formed by all nodes
    with degree k > k'. The calculation of k-density implies an iterative
    process. At each step, all nodes with degree k' < k are removed from the
    network (adjmatrix) and the density of the remaining subgraph is calculated.
    See original paper by S. Zhou and R.J. Modragon, IEEE Communication
    Letters 8(3), 180-182 (2004).
    This function returns the k-density for all values of k', from k' = 0
    (corresponding to the original graph) to k' = kmax, with kmax being the
    degree of the most connected node. It accepts weighted adjacency matrices
    but it ignores the weights of the links.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    rctype : string, optional
        Defines how to make use of the degrees depending on whether the
        network is directed or undirected.
        - 'undirected' only if the network is undirected. Raises an error if
        selected with a directed input adjmatrix.
        - 'outdegree', if the network is directed, the k-density is computed
        considering the nodes with output degree out-k' > k.
        - 'indegree', if the network is directed, the k-density is computed
        considering the nodes with input degree in-k' > k.
        - 'average', if the network is directed, the k-density is computed
        considering that the degree of the nodes is k' = 1/2 (in-k + out-k).
        Only use 'average' when the input and output degrees of the network
        are reasonably symmetric.

    Returns
    -------
    kdensity : ndarray of dtype=float
        An array of length k_max containing, for each index k, the k-density
        of the network.

    See Also
    --------
    RichClub : Identifies the subset of hubs with dense interconnectivity.
    """
    # 0) SECURITY CHECKS
    keylist = ('undirected', 'outdegree', 'indegree', 'average')
    if rctype not in keylist:
        raise KeyError("Enter a valid rctype:", keylist)

    # Convert the network in binary
    adjmatrix = adjmatrix.astype('bool')

    # Select the proper data
    indegree, outdegree = Degree(adjmatrix, True)
    if rctype == 'undirected':
        if Reciprocity(adjmatrix) < 1.0:
            raise TypeError("Option 'undirected' requires an undirected adjacency matrix")
        degree = outdegree
    elif rctype == 'outdegree':
        degree = outdegree
    elif rctype == 'indegree':
        degree = indegree
    elif rctype == 'average':
        degree = 0.5 * (indegree + outdegree)

    # 1) Prepare for calculations
    N = len(adjmatrix)

    kmax = int( np.round((degree+0.000001).max()) )
    kdensity = np.zeros(kmax+1, np.float64)

    # Density of the original network
    initialL = adjmatrix.sum()
    kdensity[0] = np.float64(initialL) / (N*(N-1))

    # 2) Compute the k-density of all degrees
    klist = np.unique( np.round((degree+0.000001)) )
    for k in range(1,kmax+1):
        # Avoid unnecessary iterations
        if k in klist:
            # 2.1) Remove the links of all nodes with degree = k
            if rctype == 'average':
                nodes = np.where(degree<=k)[0]
            else:
                nodes = np.where(degree==k)[0]
            adjmatrix[nodes] = 0
            adjmatrix[:,nodes] = 0
            degree[nodes] = 0

            # 2.2 Compute the k-density of the remainig network
            Lk = adjmatrix.sum()
            Nk = len(degree.nonzero()[0])

            if Nk > 1:
                kdensity[k] = np.float64(Lk)/(Nk*(Nk-1))

        else:
            kdensity[k] = kdensity[k-1]

    return kdensity

def RichClub(adjmatrix, kdensthreshold=0.8, rctype='undirected'):
    """Identifies the subset of hubs with dense interconnectivity.

    Conceptually, a network is said to have a rich-club when it contains hubs
    (largely connected nodes) and those hubs are densely interconnected
    forming a cluster (or community).
    The identification of a rich-club relies on the k-density, a metric that
    iteratively calculates the density of the subnetwork formed by the nodes
    with degree larger than (or equal to) k', for all k' from 0 to kmax.
    See original paper by S. Zhou and R.J. Modragon, IEEE Communication
    Letters 8(3), 180-182 (2004).

    NOTE
    Unfortunately, rich-club identification also implies a few arbitrary choices
    and no strict criteria are commonly agreed. This function returns as the
    rich-club the set of nodes remaining in the graph at the point when
    k-density first reaches a given threshold: 'kdensthreshold'. The rich-club
    will be empty if k-density does not reach the value. 'kdensthreshold' is
    an optional parameter. Initially 'kdensthreshold' is set to 0.8, which is a
    very large density, rarely achieved as the internal density by communities
    found in empirical networks by community detection methods. However, we
    recommend the user to explore different values of threshold and, specially,
    to always study the evolution of k-density as a function of degree k,
    information that is also returned by the function.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    kdensthreshold : float
        The value k-density needs to arrive, for the presence of a rich-club to
        be is considered. Real value bounded between 0 and 1.
    rctype : string, optional
        Defines how to make use of the degrees depending on whether the
        network is directed or undirected.
        - 'undirected' only if the network is undirected. Raises an error if
        selected with a directed input adjmatrix.
        - 'outdegree', if the network is directed, the k-density is computed
        considering the nodes with output degree out-k' > k.
        - 'indegree', if the network is directed, the k-density is computed
        considering the nodes with input degree in-k' > k.
        - 'average', if the network is directed, the k-density is computed
        considering that the degree of the nodes is k' = 1/2 (in-k + out-k).
        Only use 'average' when the input and output degrees of the network
        are reasonably symmetric.

    Returns
    -------
    kdensity : ndarray of dtype=float
        An array of length k_max containing, for each index k, the k-density
        of the network.
    kdecision : integer
        The degree at which k-density overcomes value 'kdensthreshold'.

    See Also
    --------
    k_Density : Iterative density of subnetworks with degree > k', from k' = 0
                to kmax.
    """
    # 0) SECURITY CHECKS
    if kdensthreshold < 0.0 or kdensthreshold > 1.0:
        raise ValueError("kdensthreshold parameter out of bounds. Please enter a value between 0 and 1.")
    keylist = ('undirected', 'outdegree', 'indegree', 'average')
    if rctype not in keylist:
        raise KeyError("Enter a valid rctype:", keylist)

    # 1) PREPARE FOR THE CALCULATIONS
    # Convert the network in binary
    adjmatrix = adjmatrix.astype('bool')

    # Select the proper data
    indegree, outdegree = Degree(adjmatrix, True)
    if rctype == 'undirected':
        if Reciprocity(adjmatrix) < 1.0:
            raise TypeError("Option 'undirected' requires an undirected adjacency matrix")
        degree = outdegree
    elif rctype == 'outdegree':
        degree = outdegree
    elif rctype == 'indegree':
        degree = indegree
    elif rctype == 'average':
        degree = 0.5 * (indegree + outdegree)

    # 1) DO THE CALCULATIONS
    N = len(adjmatrix)
    # Compute the k-density for all values of k, from 0 to kmax
    kdensity = k_Density(adjmatrix, rctype=rctype)

    if kdensity.max() >= kdensthreshold:
        # Identify the degree at which k-density overcomes the given threshold
        kdecision = np.where(kdensity >= kdensthreshold)[0][0]
        # Identify indices of the nodes forming the rich-club
        richclub = np.where(degree >= kdecision)[0]
    else:
        # If kdensity does not reach threshold for any k, return an empty set.
        kdecision = np.NaN
        richclub = np.empty(0, dtype=np.int)

    return (kdensity, kdecision, richclub)

def MatchingIndex(adjmatrix, normed=True):
    """Computes the number of common neighbours of every pair of nodes.

    The matching index of two nodes i and j is the number of common
    neighbours they are linked with.

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

    MImatrix = np.identity(N, np.float64)

    for i in range(N):
        ineighbours = set(adjmatrix[i].nonzero()[0])
        for j in range(i,N):
            jneighbours = set(adjmatrix[j].nonzero()[0])

            # Intersection and union of the sets
            mi = len(ineighbours & jneighbours)

            if normed:
                union = ineighbours | jneighbours
                norm = len(union)
                # Avoid counting the explicit links i-->j and j-->i
                if i in union: norm -= 1
                if j in union: norm -= 1
                # Normalize and save the value avoiding ZeroDivision errors
                if norm > 0:
                    mi = float(mi) / norm
                    MImatrix[i,j] = mi
                    MImatrix[j,i] = mi

            else:
                # Save the value
                MImatrix[i,j] = mi
                MImatrix[j,i] = mi

    return MImatrix

###############################################################################
"""PATHS, CYCLES AND DISTANCE FUNCTIONS"""
def FloydWarshall(adjmatrix, weighted_dist = False):
    """Computes the pathlength between all pairs of nodes in a network..

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    weighted_dist : boolean, optional
        True if the path distance shall be computed considering the weights
        of the links, False, otherwise. If the adjmatrix is a weighted
        network but'weighted = False', the unweighted graph distance is
        computed.

    Returns
    -------
    distmatrix : ndarray of rank-2
        The pairwise distance matrix dij of the shortest path between
        nodes i and j.
    """
    # Prepare for computations
    if weighted_dist:
        distmatrix = np.where(adjmatrix == 0, np.inf, adjmatrix)
    else:
        distmatrix = np.where(adjmatrix == 0, np.inf, 1)

    # Check whether the network is directed or undirected
    recip = Reciprocity(adjmatrix)

    N = len(adjmatrix)
    # Run the Floyd-Warshall algorithm - Undirected networks
    if recip == 1.0:
        for k in range(N):
            for i in range(N):
                for j in range(i,N):
                    d = distmatrix[i,k] + distmatrix[k,j]
                    if distmatrix[i,j] > d:
                        distmatrix[i,j] = d
                        distmatrix[j,i] = d

    # Run the Floyd-Warshall algorithm - directed networks
    else:
        for k in range(N):
            for i in range(N):
                for j in range(N):
                    d = distmatrix[i,k] + distmatrix[k,j]
                    if distmatrix[i,j] > d:
                        distmatrix[i,j] = d

    return distmatrix

def PathsAllinOne(adjmatrix):
    """Returns pathlength and betweenness. Finds all shortest paths and cycles.

    This function computes in a single run the following network properties:
    - The shortest pathlength between all pairs of nodes.
    - The betweennes centrality of every node in the network.
    - Finds all Hamiltonian shortest paths in the network.
    - Finds all Hamiltonian shortest cycles in the network.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.

    Returns
    -------
    distmatrix : ndarray of rank-2
        The pairwise distance matrix between all nodes.
    betweenness : ndarray
        An array containing the betweenness centrality of every node in the
        network. Result is in absolute values, i.e., the total number of
        paths in which the node participates as intermediate node.
    allpaths : dictionary
        A dictionary containing all Hamiltonian shortest paths in the
        network sorted by length. For example, allpaths[4] is a list of
        all shortest paths of length 3 in the network. Each path is
        represented as an ordered list of nodes, e.g., [0,4,3,,21,10] is a
        path from node 0 to node 10.
    allcycles : dictionary
        A dictionary containing all Hamiltonian shortest cycles in the
        network sorted by length. For example, allcycles[3] is a list of
        all shortest cycles of length 3 in the network. Each cycle is
        represented as an ordered list of nodes, e.g., [0,4,3,0] is a cycle
        of length 3.

    See Also
    --------
    FloydWarshall : Computes the pathlength between all pairs of nodes.
    AllShortestPaths : Find all shortest paths between two nodes.
    CleanCycles : Finds repeated cycles and leaves only one copy.
    CleanPaths : Removes opposite running paths from undirected graphs.

    Warnings
    --------
    1. If the adjmatrix is undirected every cycle is found and saved n times,
    where n is the number of nodes in the cycle, e.g., [1,3,4], [4,1,3] and
    [3,1,4]. If adjmatrix is directed, cycles might appear an unpredictable
    number of times or only once. Use function CleanCycles() in
    galib.tools module to remove repeated cycles.
    2. If adjmatrix is undirected every path is saved twice, in the opposite
    direction, e.g. [1,3,4] and [4,3,1]. Use CleanPaths() in galib.tools module
    to remove duplicated paths from undirected graphs.
    """
    # 0) PREPARE FOR THE CALCULATIONS
    N = len(adjmatrix)
    distmatrix = np.where(adjmatrix, 1, np.inf)
    betweenness = np.zeros(N, np.int)
    allpaths = {}
    allcycles = {}

    # 1) FIND PATHS OF LENGTH 1, CONVERT NETWORK INTO DICTIONARY
    # The network is converted into dictionary to improve speed of step 2.1)
    # at the cost of the memory it costs to host the network again.
    allpaths[1] = []
    dicnet = {}
    for i in range(N):
        neighs = adjmatrix[i].nonzero()[0]
        dicnet[i] = neighs.tolist()
        for j in neighs:
            allpaths[1].append([i,j])
            distmatrix[i,j] = 1

    # 2) FIND THE REST OF THE PATHS
    for length in range(2,N):

        # 2.0) Check if the algorithm has already finished. If so, end
        if not allpaths[length-1]:
            del allpaths[length-1]
            break

        allpaths[length] = []
        # 2.1) To all known paths of length - 1, try to append a new node
        for path in allpaths[length-1]:
            start = path[0]
            end = path[-1]
            for j in dicnet[end]:

                # 2.2) If new path + [j] is a cycle, treat separately
                if j == start:
                    idx = np.argmin(path)
                    cyc = path[idx:] + path[:idx]
                    if length in allcycles:
                        if cyc in allcycles[length]:
                            continue
                        else:
                            allcycles[length].append(cyc)
                    else:
                        allcycles[length] = [cyc]

                    # Update the distance matrix, only if necessary
                    if distmatrix[start,start] >= length:
                        distmatrix[start,start] = length
                    continue

                # 2.3) Discard path if it is longer than known dist(start,j)
                if distmatrix[start,j] < length: continue

                # 2.4) Discard non-Hamiltonian paths
                if j in path: continue

                # 2.5) Only if path+[j] is a new Hamiltonian path, include it
                allpaths[length].append(path+[j])
                # Update the distance matrix
                distmatrix[start,j] = length
                # Update the betweenness centrality
                for node in path[1:]:
                    betweenness[node] += 1

    # 3) CORRECT FOR UNDIRECTED GRAPHS
    # For the moment I only correct the betweenness, but I should remove
    # also all undirected paths that run inversely and are duplicated.
    recip = Reciprocity(adjmatrix)
    if recip == 1:
        betweenness = (0.5*betweenness).astype(np.int)


    # 4) CLEAN TRASH AND FINISH
    del dicnet
    return distmatrix, betweenness, allpaths, allcycles

def ShortestPaths(adjmatrix, start, end, length, queue = [], paths = []):
    """Finds all the shortest paths between two nodes.

    A recursive function that uses the deep-first-search (DFS) strategy to
    navigate the graph starting at node 'start' and finishing the search at
    branches of length 'length'. adjmatrix can be a weighted adjacency
    matrix but weights of links will be ignored.

    Usage
    -----
    Always call the function with 'queue' being the start node AND with
    parameter 'paths=[]' typed. This is very important due to the recursive
    nature of the algorithm. If 'paths' is left empty, the paths obtained in
    two consecutive calls will accumulate. Proper calling example:

    result = AllHamiltonianPaths(network,node1,node2,queue=[node1],paths=[])

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    start : integer
        The node from which the tree search starts.
    end : integer
        The target node at which the paths finish.
    length : integer
        The desired graph distance between nodes 'start' and 'end'. To find
        all shortest paths between them, 'length' has to be the shortest graph
        distance which must have been computed beforehand by other methods,
        e.g. Dijkstra's or the Floyd-Warshall algorithm.
    queue : list of integers
        The current path of nodes that the function explores. SEE 'Usage'!!!

    Returns
    -------
    paths : lists of lists
        The list of all the shortest paths between nodes 'start' and 'end'.
        LEAVE EMPTY at function call!! See Usage.

    See Also
    --------
    FloydWarshall : Computes all-to-all shortest graph distance
    """
    # 1) Enqueue a neighbour of the last node in queue
    # Only if it is not already in queue (Hamiltonian paths)
    for j in adjmatrix[start].nonzero()[0]:
        if j in queue[1:]: continue
        if len(queue) > length: continue

        queue.append(j)
        ShortestPaths(adjmatrix, j, end, length, queue, paths)

        # 2) Check if current queue is a path of desired characteristics
        if (len(queue)-1 == length and j == end):
            paths.append(list(queue))

        # 3) Remove last item in queue
        queue.pop()

    return paths


############################################################################
"""COMPONENTS, COMMUNITIES, K-CORES..."""
def ConnectedComponents(distmatrix, directed=False, showall=True):
    """Finds all the connected components in a network out of a distance
    matrix.

    A strongly connected component is a set of nodes for which there is
    at least one path connecting every two nodes within the set.
    The function works both for directed and undirected networks, provided
    the adequate distance matrix is given.

    Parameters
    ----------
    distmatrix : ndarray of rank-2
        The pairwise graph distance matrix of the network, usually the
        output of function FloydWarshall().
    directed : boolean, optional
        'True' if the network is directed, 'False' if it is undirected.
    showall : boolean, optional
        If 'True' the function returns all strong components, including
        independent nodes. If 'False' it returns only components of two
        or more nodes.

    Returns
    -------
    components : list
        A list containing the components as ordered lists of nodes.

    See Also
    --------
    FloydWarshall : Pairwise graph distance between all nodes of a network.
    """
    N = len(distmatrix)

    # 1) Detect nodes that are connected in both directions
    newmatrix = np.where(distmatrix < N, 1, 0)

    # If network is directed, consider only pairs with a reciprocal path
    if directed:
        newmatrix = newmatrix * newmatrix.T

    # 2) Sort the nodes into their components
    nodelist = range(N)
    components = []
    while nodelist:
        # Take the first node. This helps keeping the output sorted
        node = nodelist[0]
        if newmatrix[node,node]:
            # Find the component to which the node belongs to
            comp = list(newmatrix[node].nonzero()[0])
            components.append(comp)
            # Remove nodes in comp from nodelist
            for neigh in comp:
                nodelist.remove(neigh)
            # Clean trash
            del comp
        else:
            # The node is independent. Remove from list and continue
            if showall:
                components.append([node])
            nodelist.remove(node)
            continue

    del newmatrix
    return components

def AssortativityMatrix(adjmatrix, partition, norm=None, maxweight=1.0):
    """Returns the assortativity matrix of network given a partition of nodes.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays). 'partition' may contain any arbitrary grouping of network
        nodes, e.g., overlapping subsets are also accepted.
    norm : boolean or string, optional
        Defines the kind of normalization for the counts of links:
        - If norm=None, returns the number of links between subsets of nodes.
        - If norm='linkfraction', returns the number of links divided by the
        total number of links in the network.
        - If norm='linkprobability' returns the number of links between two
        subsets, divided by the total number of possible links between them.
    maxweight : floating-point scalar, optional
        Largest possible weight of the links.

    Returns
    -------
    assortmatrix : ndarray of rank-2 and dtype=float
        Assortativity matrix of shape Nc x Nc, where Nc is the number of
        subsets of nodes in 'partition'.

    Notes
    -----
    The function accepts weighted adjacency matrices but assumes link
    weights to lie between 0 and 'maxweight'. See documentation.

    See Also
    --------
    ParticipationMatrix : Probability of nodes to belong to a community.
    """
    # Security check
    keylist = [None, 'linkfraction', 'linkprobability']
    if norm not in keylist:
        raise KeyError("Enter a valid norm:", keylist)

    N = len(adjmatrix)
    Ncoms = len(partition)

    # Calculate the assortativity matrix
    assortmatrix = np.zeros((Ncoms,Ncoms), np.float)

    if norm == 'linkprobability':
        for c1 in range(Ncoms):
            com1 = partition[c1]
            for c2 in range(Ncoms):
                com2 = partition[c2]
                submat = tools.ExtractSubmatrix(adjmatrix, com1, com2)
                assortmatrix[c1,c2] = submat.sum()
                # Normalise, avoiding self-loops
                if c1 == c2:
                    ncom1 = len(com1)
                    assortmatrix[c1,c2] /= (maxweight * ncom1*(ncom1-1))
                else:
                    assortmatrix[c1,c2] /= (maxweight * len(com1) * len(com2))

    else:
        for c1 in range(Ncoms):
            for c2 in range(Ncoms):
                submat = tools.ExtractSubmatrix(adjmatrix, partition[c1], \
                                                  partition[c2])
                assortmatrix[c1,c2] = submat.sum()

        if norm == 'linkfraction' and maxweight == 1.0:
            assortmatrix /= adjmatrix.sum()
        elif norm == 'linkfraction' and maxweight != 1.0:
            L = len(adjmatrix.flatten().nonzero()[0])
            assortmatrix /= (maxweight*L)

    return assortmatrix

def Modularity(adjmatrix, partition):
    """Computes the Newman modularity given a partition of nodes.

    It computes modularity for both weighted or unweighted and for directed
    or undirected networks after a generalization of the modularity measure
    by S. Gomez, P. Jensen & A. Arenas [Phys. Rev. E 80,016114 (2009)].
    See Notes and the GAlib documentation.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    Q : float scalar
        The modularity value of the network for the given partition

    Notes
    -----
    If 'adjmatrix' is the binary adjacency matrix of the network and 'degree'
    the usual degree, the function returns the unweighted modularity.
    If 'adjmatrix' is a weighted adjacency matrix and 'degree' its weighted
    degree or intensity, it returns the weighted modularity.
    The algorithm authomatically computes modularity for both directed and
    undirected networks according to input 'matrix'.
    """
    # Prepare for calculations
    N = len(adjmatrix)
    Ncoms = len(partition)
    L = adjmatrix.sum()

    indegree, outdegree = Intensity(adjmatrix, directed=True)

    # Compute the modularity
    Q = 0.0
    L_norm = 1./L
    for s, community in enumerate(partition):
        submat = tools.ExtractSubmatrix(adjmatrix, community)
        # Add the fraction of internal links
        Q += float(submat.sum()) * L_norm
        # Minus the expected fraction of links
        productsubmat = np.outer(outdegree[community], indegree[community])
        Q -= productsubmat.sum() * L_norm**2

    return Q

def K_Core(adjmatrix, kmin):
    """Finds the K-core of a network with degree k >= kmin.

    A k-core is the largest subgraph for which all nodes have degree >= k
    within the subgraph. Hence, a k-core is composed of at least k+1 nodes.

    Parameters
    ----------
    adjmatrix : ndarray
        The adjacency matrix of the network.
    kmin : integer
        The degree for which the k-core is desired.

    Returns
    -------
    nodelist : list of integers
        The k-core of the network where k = kmin.

    Usage
    -----
    The function accepts weighted networks, but ignores the link weights.
    In case of directed networks, for the k-core decomposition of output
    degrees pass the usual adjacency matrix such that Aij = 1 if there
    is a link from i to j. For the decomposition based on input degree,
    pass the transposed adjacency matrix to the function.

    See Also
    --------
    K_Shells : Returns the k-shells of a network for all k.
    """
    # Prepare for calculations
    adjmatrix = adjmatrix.copy()
    N = len(adjmatrix)
    nodelist = set(np.arange(N))
    degree = Degree(adjmatrix)

    # Start detecting the kmin-core
    done = False
    while not done:
        # 1) Select nodes with degree <= kmin and include them in the shell
        nodes = np.where(degree < kmin)[0]
        nodes = list(set(nodes) & nodelist)

        # 2) Remove the selected nodes from the network and from 'nodelist'
        adjmatrix[nodes] = 0
        adjmatrix[:,nodes] = 0
        for node in nodes:
            nodelist.remove(node)

        # 3) Find if the removal caused other nodes to have degree <= kmin.
        # In last iteration NonZeromin() raises a ValueError if
        # 'nodelist' is already empty before end of loop.
        degree = Degree(adjmatrix)
        try:
            newkmin = tools.NonZeroMin(degree)
        except ValueError:
            done = True

        # If all remaining nodes have degree > kmin, finish,
        # otherwise, continue updating the core.
        if newkmin >= kmin:
            done = True

    return list(nodelist)

def K_Shells(adjmatrix):
    """Returns the K-shells of a network for all k from kmin to kmax.

    A k-core is the largest subgraph for which all nodes have degree >= k
    within the subgraph. Hence, a k-core is composed of at least k+1 nodes.
    A k-shell is composed by all nodes with coreness k, such that the
    k-core is the union of all k'-shells with k' >= k

    Parameters
    ----------
    adjmatrix : ndarray
        The adjacency matrix of the network.

    Returns
    -------
    kshells : dictionary
        Dictionary containing the list of nodes forming k-shell for
        each key k.

    Usage
    -----
    The function accepts weighted networks, but ignores the link weights.
    In case of directed networks, for the k-core decomposition of output
    degrees pass the usual adjacency matrix such that Aij = 1 if there
    is a link from i to j. For the decomposition based on input degree,
    pass the transposed adjacency matrix to the function.

    See Also
    --------
    K_Core : Returns the core of a network for a given degree threshold k.
    """
    # 0) Prepare for calculations
    adjmatrix = adjmatrix.copy()
    N = len(adjmatrix)
    nodelist = set(np.arange(N))

    # Find the smallest non-zero degree to start from
    degree = Degree(adjmatrix)
    kmin = tools.NonZeroMin(degree)

    # 1) Start computing the k-shells
    kshells = {}
    while nodelist:
        shell = []
        done = False
        # Find the next shell
        while not done:
            # 1.1) Select nodes with degree <= kmin and include them in shell
            nodes = np.where(degree<=kmin)[0]
            nodes = list(set(nodes) & nodelist)
            shell += list(nodes)

            # 1.2) Remove the selected nodes from network and from 'nodelist'
            adjmatrix[nodes] = 0
            adjmatrix[:,nodes] = 0
            for node in nodes:
                nodelist.remove(node)

            # 1.3) Find if removal caused other nodes to have degree <= kmin
            # In last iteration NonZeromin() raises a ValueError if
            # 'nodelist' is already empty before end of loop
            degree = Degree(adjmatrix)
            try:
                newkmin = tools.NonZeroMin(degree)
            except ValueError:
                kshells[kmin] = shell
                done = True

            # If all remaining nodes have degree > kmin, look for next shell
            # otherwise, continue updating the current shell
            if newkmin > kmin:
                kshells[kmin] = shell
                kmin = newkmin
                done = True

    return kshells


######################################################################
"""ROLES OF NODES IN NETWORKS WITH COMMUNITY (ASSORTATIVE) ORGANIZATION"""
def GlobalHubness(adjmatrix):
    """Computes the global hubness of all nodes in a network.

    Hubness is the degree of a node weighted by the expected degree
    distribution in random graphs of same size and density. See Equation (4)
    of Klimm et al. New J. Phys. 16:125006 (2014).

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.

    Returns
    -------
    globalhubness : ndarray (float64)
        Global hubness of every node.

    See Also
    --------
    LocalHubness : Given a partition, computes the local hubness of all nodes.
    NodeParticipation : Participation index of every node given a partition of the network.
    NodeDispersion : Dispersion index of every node given a partition of the network.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    ParticipationVectors : Computes the probability of nodes to belong to every community.
    NodeRoles :

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)

    dens = Density(adjmatrix)
    invnorm = 1. / np.sqrt((N-1) * dens * (1.0 - dens))
    degree = Degree(adjmatrix)

    globalhubness = invnorm * (degree - (N-1)*dens)
    return globalhubness

def LocalHubness(adjmatrix, partition):
    """Given a partition, computes the local hubness of all nodes.

    Hubness is the degree of a node weighted by the expected degree
    distribution in random graphs of same size and density. e Equation (4)
    of Klim et al. New J. Phys. 16:125006 (2014).
    Local hubness is the hubness applied to the subgraph formed by the
    community the node belongs to.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    localhubness : ndarray (float64)
        Local hubness of every node.

    See Also
    --------
    GlobalHubness : Hubness of nodes within their community.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    ParticipationVectors : Computes the probability of nodes to belong to every community.
    NodeParticipation : Participation index of every node given a partition of the network.
    NodeDispersion : Dispersion index of every node given a partition of the network.
    NodeRoles :

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)
    ncoms = len(partition)

    localhubness = np.zeros(N,np.float64)
    for n, com in enumerate(partition):

        # Skip nodes of only one node
        if len(com) == 1: continue

        # Compute the hubness of nodes in the isolated community
        subnet = tools.ExtractSubmatrix(adjmatrix,com)
        hubness = GlobalHubness(subnet)

        if np.isnan(hubness.min()): continue

        localhubness[com] = hubness

    return localhubness

def ParticipationMatrix(adjmatrix, partition):
    """Given a partition of the network, it returns the participation matrix.

    A matrix of shape N x n, where N is the number of nodes and n is the
    number of communities. Elements a(i,s) of the matrix are the number of
    neighbours (internal degree) that node i has in community s.

    Parameters
    ----------
    adjmatrix : ndarray
        The adjacency matrix of the network.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    pmatrix : ndarray of rank-2 and shape N x n.

    See Also
    --------
    GlobalHubness : Computes the global hubness of all nodes in a network.
    LocalHubness : Given a partition, computes the local hubness of all nodes.
    PArticipationVectors : Computes the probability of nodes to belong to every community.
    NodeParticipation : Participation index of every node given a partition of the network.
    NodeDispersion : Dispersion index of every node given a partition of the network.
    RolesNode :

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)
    ncomms = len(partition)
    partitionmatrix = np.zeros((N,ncomms), np.uint)

    # 1) CONSTRUCT THE PARTITION MATRIX, S (1 if node in module c, 0 otherwise)
    for c in range(ncomms):
        partitionmatrix[partition[c],c] = 1

    # 2) COMPUTE THE PARTICIPATION MATRIX
    # adjmatrix.astype(bool) for cases in which adjmatrix is weighted
    pmatrix = np.dot(adjmatrix.astype(bool), partitionmatrix)

    return pmatrix

def ParticipationVectors(adjmatrix, partition):
    """Computes the probability of nodes to belong to every community.

    The probability of node i to belong to community c is estimated as the
    fraction of nodes in c with which i is connected to: k_ic / N_c, where
    k_ic is the degree of i in c and N_c is the size of the community.
    If i belongs to community c, then the norm is (N_c -1). Finally, the
    fractions are normalised such that their sum is 1.

    Parameters
    ----------
    adjmatrix : ndarray
        The adjacency matrix of the network.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    pmatrix : ndarray of rank-2 and shape N x n.
        Every row contains the likelihood of the node to belong to each
        of the communities.

    See Also
    --------
    GlobalHubness : Computes the global hubness of all nodes in a network.
    LocalHubness : Given a partition, computes the local hubness of all nodes.
    NodeParticipation : Participation index of every node given a partition of the network.
    NodeDispersion : Dispersion index of every node given a partition of the network.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    RolesNode :

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)
    ncomms = len(partition)
    commsizes = np.zeros(ncomms, np.float64)
    partitionmatrix = np.zeros((N,ncomms), np.uint8)

    # 1) COMPUTE FIRST THE PARTICIPATION MATRIX
    # 1.1) Construct the partition matrix
    for c in range(ncomms):
        partitionmatrix[partition[c],c] = 1
        commsizes[c] = len(partition[c])

    # 1.2) Compute the participation matrix
    # adjmatrix.astype(bool) for cases in which adjmatrix is weighted
    pmatrix = np.dot(adjmatrix.astype(bool), partitionmatrix)

    # 2) NOW COMPUTE THE PARTICIPATION VECTORS
    pmatrix = pmatrix.astype(np.float64)
    for i in range(N):
        # 2.1) The fraction of neighbours node i connects with,
        # in every module
        for c in range(ncomms):
            if partitionmatrix[i,c]:
                pmatrix[i,c] /= (commsizes[c] - 1.0)
            else:
                pmatrix[i,c] /= commsizes[c]

        # 2.2) Finally, normalize the vector to sum 1
        pmatrix[i] /= pmatrix[i].sum()

    return pmatrix

def NodeParticipation(adjmatrix, partition):
    """Participation index of every node given a partition of the network.

    Computes the participation index of every node. If a node is only
    connected to other nodes of the same community, then pi = 0. If the
    node is equally likely connected with ALL the communities, then pi = 1.
    For a detailed definition see Equation (7) of Klimm et al. New J. Phys.
    16:125006 (2014).

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    nodeparticip : ndarray of dtype=Float64 and rank-N
        The participation index of each node.

    See Also
    --------
    GlobalHubness : Computes the global hubness of all nodes in a network.
    LocaHubness : Given a partition, computes the local hubness of all nodes.
    NodeDispersion : Dispersion index of every node given a partition of the network.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    ParticipationVectors : Computes the probability of nodes to belong to every community.

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    # 1) Compute first the participation vectors
    particvectors = ParticipationVectors(adjmatrix,partition)

    # 2) Reduce the vector of each node into a single normalised scalar
    ncomms = np.float64(len(partition))
    nodeparticip = 1. - ncomms / np.sqrt(ncomms - 1) * particvectors.std(axis=1)

    return nodeparticip

def NodeDispersion(adjmatrix, partition):
    """Dispersion index of every node given a partition of the network.

    Computes the dispersion index of every node. Dispersion is a measure to
    quantify whether a node is well classified within the partition, or it
    is missaclassified. If a node is only connected to nodes in the same
    community, then di = 0. If the node is equally likely connected with two
    of more communities, then pi = 1.
    For a detailed definition see Equation (7) of Klimm et al. New J. Phys.
    16:125006 (2014).

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    nodedispersion : ndarray of dtype=Float64 and rank-N
        The participation index of each node.

    See Also
    --------
    GlobalHubness : Computes the global hubness of all nodes in a network.
    LocalHubness : Given a partition, computes the local hubness of all nodes.
    NodeParticipation : Participation index of every node given a partition of the network.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    ParticipationVectors : Computes the probability of nodes to belong to every community.

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)
    # 1) Compute first the participation vectors
    particvectors = ParticipationVectors(adjmatrix,partition)

    # 2) Reduce the vector of each node into a single normalised scalar
    nodedispersion = np.zeros(N, np.float64)
    for i in range(N):
        idx = particvectors[i].nonzero()[0]
        vector = particvectors[i,idx]
        ncomms = np.float64(len(vector))
        if ncomms > 1:
            nodedispersion[i] = 1. - ncomms / np.sqrt(ncomms - 1) * vector.std()

    return nodedispersion

def RolesNodes(adjmatrix, partition):
    """Computes all four parameters to characterise the roles of nodes.

    Following the definitions in Klimm et al. New J. Phys. 16:125006 (2014),
    it computes the four parameters that characterise the roles that nodes
    take in a network given a partition of its nodes into communites, classes
    or any predefined categories.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network. Weighted links are ignored.
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    globalhubness : ndarray (float64)
        Global hubness of every node.
    localhubness : ndarray (float64)
        Local hubness of every node.
    nodeparticip : ndarray of dtype=Float64 and rank-N
        The participation index of each node.
    nodedispersion : ndarray of dtype=Float64 and rank-N
        The participation index of each node.

    See Also
    --------
    GlobalHubness : Computes the global hubness of all nodes in a network.
    LocalHubness : Given a partition, computes the local hubness of all nodes.
    NodeParticipation : Participation index of every node given a partition of the network.
    ParticipationMatrix : Given a partition of the network, it returns the participation matrix.
    ParticipationVectors : Computes the probability of nodes to belong to every community.

    Citation
    --------
    F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths & G. Zamora-Lopez,
    "Individual nodeʼs contribution to the mesoscale of complex networks."
    New Journal of Physics 16:125006 (2014).

    """
    N = len(adjmatrix)
    ncomms = len(partition)

    # 1) COMPUTE THE GLOBAL HUBNESS
    dens = Density(adjmatrix)
    invnorm = 1. / np.sqrt((N-1) * dens * (1.0 - dens))
    degree = Degree(adjmatrix)

    globalhubness = invnorm * (degree - (N-1)*dens)

    # 2) COMPUTE THE LOCAL HUBNESS AND CREATE THE PARTITION MATRIX
    localhubness = np.zeros(N,np.float64)
    commsizes = np.zeros(ncomms,np.float64)
    partitionmatrix = np.zeros((N,ncomms), np.uint)
    for c, com in enumerate(partition):
        commsizes[c] = len(partition[c])
        partitionmatrix[com,c] = 1

        # Skip nodes of only one node
        if len(com) == 1: continue

        # Compute the hubness of nodes in the isolated community
        subnet = tools.ExtractSubmatrix(adjmatrix,com)
        hubness = GlobalHubness(subnet)

        if np.isnan(hubness.min()): continue
        localhubness[com] = hubness

    # 3) COMPUTE THE PARTICIPATION VECTORS
    particvectors = np.dot(adjmatrix.astype(bool), partitionmatrix)
    particvectors = particvectors.astype(np.float64)
    for i in range(N):
        for c in range(ncomms):
            if partitionmatrix[i,c]:
                particvectors[i,c] /= (commsizes[c] - 1.0)
            else:
                particvectors[i,c] /= commsizes[c]

        # Finally, normalize the vector to sum 1
        particvectors[i] /= particvectors[i].sum()

    # 4) COMPUTE THE PARTICIPATION AND DISPERSION INDICES
    # The participation index
    norm = np.float64(ncomms) / np.sqrt(ncomms - 1.0)
    nodeparticip = 1. - norm * particvectors.std(axis=1)

    # The dispersion index
    nodedispersion = np.zeros(N, np.float64)
    for i in range(N):
        idx = particvectors[i].nonzero()[0]
        vector = particvectors[i,idx]
        ncomms = np.float64(len(vector))
        if ncomms > 1:
            nodedispersion[i] = 1. - ncomms / np.sqrt(ncomms - 1) * vector.std()

    return globalhubness, localhubness, nodeparticip, nodedispersion

def ParticipationIndex_GA(participmatrix):
    """Returns the participation index as defined by Guimera & Amaral.

    Given a partition of the network into communities, the participation
    index quantifies how much are the links of a node distributed along
    all the communities of the network. This function computes the definition
    given by Guimera & Amaral, J. Stat. Mech. P02001 (2011).

    Parameters
    ----------
    participmatrix : ndarray
        A matrix of shape N x n, where N is the number of nodes and n is the
        number of communities. Elements, a_is, of the matrix are the number of
        neighbours (degree) that node i has in community s.

    Returns
    -------
    participindex : ndarray of rank-1 and size N
        Participation index of every node.

    Notes
    -----
    The participation index, as originally defined in Guimera & Amaral,
    J. Stat. Mech. P02001 (2011), is misleading. Use only for comparative
    reasons. The index is defined to be between 0 and 1, taking value 0 when
    a node has only connections in its own community (peripheral node),
    and 1 when its links are equivalently distributed among, so the node is
    unclassifiable in the partition (kinless). However, the true range of this
    participation index depends on the number of communities in the partition.

    See Also
    --------
    ParticipationMatrix : The number of neighbours of a node in all communities.
    ParticipationIndex : Returns the participation index of all nodes.
    LocalHubness_GA : Returns the z-score of node's local degree.
    """
    N, npart = np.shape(participmatrix)
    participindex = np.zeros(N, np.float)

    for i in range(N):
        a = participmatrix[i].sum()
        if a > 0: norm_i = 1./a
        else: norm_i = 0.
        for c in range(npart):
            participindex[i] += (norm_i*participmatrix[i][c])**2

    participindex = 1. - participindex
    return participindex

def Hubness_GA(participmatrix, partition):
    """Returns the within-module degree defined by Guimera & Amaral.

    The within-module degree is a measure of local hubness. Given a network
    and a partition of it nodes into communities, the within-module degree
    is the z-score of the number of neighbours that a node has, compared with
    the degree of the other nodes in the community. See Guimera & Amaral,
    J. Stat. Mech. P02001 (2011).

    Parameters
    ----------
    participmatrix : ndarray of rank-2
        Elements a(i,s) of the matrix are the number of neighbours (internal
        degree) that node i has in community s. Is the output of function
        ParticipationMatrix().
    partition : list, tuple or array_like
        A sequence of subsets of nodes given as sequences (lists, tuples or
        arrays).

    Returns
    -------
    zscore : ndarray of rank-1
        The within-module degree of every node in the network.

    See Also
    --------
    ParticipationMatrix : The number of neighbours of a node in all communities.
    ParticipationIndex_GA : Returns the participation index (based on Guimera
        & Amaral's definiton) of all nodes given a partition.
    ParticipationIndex : Returns the participation index of all nodes given
        a partition.
    """
    N, ncoms = np.shape(participmatrix)

    zscore = np.zeros(N, np.float)
    for s, community in enumerate(partition):
        klist = participmatrix[community,s]
        avklist = klist.mean()
        devklist = klist.std()

        zscore[community] = (participmatrix[community,s] - avklist) / devklist

    return zscore

#
