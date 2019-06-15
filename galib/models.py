# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
SYNTHETIC NETWORK GENERATORS
============================

This module contains functions to generate typical synthetic networks,
including random networks and methods to rewire networks.

RANDOM NETWORK GENERATORS
-------------------------
Lattice1D
    Generates regular ring lattices (all nodes have same degree).
Lattice1D_FixLinks
    Generates ring lattices with desired number of links.
WattsStrogatzGraph
    Generates small-world networks as in the Watts & Strogatz model.
ErdosRenyiGraph
    Generates random graphs following the Erdos & Renyi model.
RandomGraph
    Generates random graphs with N nodes and L links.
BarabasiAlbertGraph
    Generates scale-free networks after the Barabasi & Albert model.
ScaleFreeGraph
    Generates scale-free graphs of given size and exponent.

NETWORK REWIRING/RANDOMIZATION ALGORITHMS
-----------------------------------------
RewireNetwork
    Randomises an input graph conserving the degrees of its nodes.
ModularityPreservingGraph
    Randomises an input graph conserving its modular structure.

HIERARCHICAL AND MODULAR (HM) NETWORK MODELS
--------------------------------------------
ModularHeterogeneousGraph
    Generates random modular networks of given module sizes and densities.
HMpartition
    Returns a partition of nodes for a hierarchical and modular network.
HMRandomGraph
    Generates random hierarchical and modular networks of desired number
    of hierarchical levels and modules.
HMCentralisedGraph
    Generates random hierarchical and modular networks of desired number
    of hierarchical levels and modules, with centralised inter-modular
    connectivity through local hubs.
RavaszBarabasiGraph
    Generates hierarchical networks following the Ravasz & Barabasi model.


...moduleauthor:: Gorka Zamora-Lopez <galib@zamora-lopez.xyz>

"""
from __future__ import division, print_function, absolute_import

import types
import numpy as np
import numpy.random
from .metrics import Reciprocity
from .tools import ExtractSubmatrix


############################################################################
"""RANDOM NETWORK GENERATORS"""
def Lattice1D(N,z):
    """Generates regular ring lattices.

    Each node is connected to its 2z closest neighbours (z on the left and z on
    the right). Hence, all nodes have same degree.

    Parameters
    ----------
    N : integer
        Size of the network (number of nodes).
    z : integer
        Diameter of the neighbourhood that every node connects with.
        For a given z, every node has degree 2*z in the resulting lattice.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the 1-dimensional lattice.

    See Also
    --------
    Lattice1D_FixLinks : Generates ring lattices with desired number of links.
    """
    # 0) SECURITY CHECK
    if z > N//2:
        raise ValueError("Largest possible z = N/2 =", N//2 )

    if z == 0:
        return np.zeros((N,N), np.uint8)

    # 1) CREATE THE LATTICE
    adjmatrix = np.zeros((N,N), np.uint8)

    # 1.1) Create the first row according to the number of neighbours
    adjmatrix[0,1:z+1] = 1
    adjmatrix[0,-z:] = 1

    # 1.2) Use numpy.roll() to copy rotated version of the first row
    for i in range(1,N):
        adjmatrix[i] = np.roll(adjmatrix[0],i)

    return adjmatrix

def Lattice1D_FixLinks(N,L):
    """Generates a 1D ring lattice with given number of links.

    Because the total number of links L is fixed, the resulting ring is not a
    perfect regular lattice in which all nodes have exactly the same degree.
    The result is quasi-regular. The largest degree found will be kmax = <k> + 2
    and the smallest degree kmin = <k> - 2.

    Parameters
    ----------
    N : integer
        Size of the network (number of nodes).
    L : integer
        Number of links the resulting network must have.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the 1-dimensional lattice.

    See Also
    --------
    Lattice1D : Generates regular ring lattices (all nodes have same degree)
    """
    # 0) SECURITY CHECK
    if L > int(0.5*N*(N-1)):
        raise ValueError( "Largest number of links 1/2 * N*(N-1) =", int(0.5*N*(N-1)) )

    if L == 0:
        return np.zeros((N,N), np.uint8)

    adjmatrix = np.zeros((N,N), np.uint8)

    # 1.2) Use numpy.roll() to copy rotated version of a pattern row
    counter = 0
    finished = False
    for k in range(1,N):
        if finished: break
        row = np.zeros(N,np.uint8)
        row[k] = 1

        for i in range(N):
            adjmatrix[i] += np.roll(row,i)
            counter += 1
            if counter == L:
                finished = True
                break

    return adjmatrix + adjmatrix.T

def WattsStrogatzGraph(N, z, prew, lattice=None):
    """Generates small-world networks as in the Watts & Strogatz model.

    See Watts & Strogatz, Nature 393(4) (1998) for details of the model.

    The algorithm rewires, clockwise, the links of the first neighbours,
    of the second neighbours, etc. with probability prew. If a link is
    chosen for rewiring, it can't create a self-loop nor a multiple edge.

    Parameters
    ----------
    N : integer
        Size of the network (number of nodes).
    z : integer
        Diameter of the neighbourhood that every node connects with.
        For a given z, every node has degree k = 2*z in the resulting lattice.
    prew : float, between 0 and 1.
        Probability that links of the 1D lattice to be rewired.
    lattice : ndarray of rank-2, optional.
        The adjacency matrix of a 1D lattice that has to be rewired.
        When several realizations of the model are desired, the initial
        lattice, usually the output of Lattice1D() function, can be passed
        what avoids having to create a new lattice for every realization and
        saving computational time. Recommended

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the rewired 1-dimensional lattice.

    Notes
    -----
    1) Unlike the Watts & Strogatz model in which k is the total number of
    neighbours of a node, here z represents the largest z-neighbourhood
    with which a node is connected (left and right). Hence, comparing the
    two parameters, k = 2*z.

    2) When several realizations of the model are desired, the initial
    lattice, usually the output of Lattice1D() function, can be passed
    to the optional argument 'lattice'. This avoids having to create a
    new lattice for every realization and saves computational time.
    Recommended for several realizations of large networks.

    3) The present algorithm is only valid for graphs. In order to obtain
    directed small-world networks, use the RewireNetwork() function on a
    previously created ring-lattice.

    See Also
    --------
    Lattice1D : Generates a 1D lattice ring network.
    RewiredNetwork : Randomized the links of a network with equal probability.
    """

    # 0) SECURITY CHECKS
    if (prew < 0 or prew > 1):
        raise ValueError( "Probability 'prew' out of bounds. Insert value between 0 and 1" )

    # 1) CREATE OR CHECK THE VALIDITY OF THE INITIAL RING-LATTICE
    # if type(lattice) == types.NoneType:
    if lattice == None:
        if z > N//2:
            raise ValueError("Largest possible z = N/2 =", N//2 )
        adjmatrix = Lattice1D(N,z)
    else:
        Nlatt = len(lattice)
        if N != Nlatt:
            raise ValueError( "N and size of given lattice od not match." )
        adjmatrix = lattice.copy()

    # 2) REWIRE THE LINKS CLOCKWISE AND BY RANK OF NEIGHBOURHOOD
    # For each rank of neighbourhood...
    for k in range(1,z+1):
        # Choose every node in clockwise direction
        for i in range(N):
            # 2.1) Rewire the link with neighbour i+k with probability prew
            if numpy.random.rand() <= prew:
                done = False
                while not done:
                    target = int( N * numpy.random.rand() )
                    # Avoid self-loops and double links
                    if target == i: continue
                    if adjmatrix[i,target]: continue

                    # Else, remove the old link and place the new one.
                    if i+k < N:
                        adjmatrix[i,i+k] = 0
                        adjmatrix[i+k,i] = 0
                    else:
                        adjmatrix[i,i+k-N] = 0
                        adjmatrix[i+k-N,i] = 0

                    adjmatrix[i,target] = 1
                    adjmatrix[target,i] = 1
                    done = True

    return adjmatrix

def ErdosRenyiGraph(N, p, directed=False, selfloops=False, outdtype=np.uint8):
    """Generates random graphs following the Erdos & Renyi model.

    In an Erdos-Renyi graph a link happens with probability p. Therefore,
    different relizations of the graph may have different number of links.
    In the ensemble average, the expected number of links of ER graphs
    is <L> = 1/2 * p * N**2 if self-loops are permitted, or
    <L> = 1/2 * p * N * (N-1) if no self-loops are permitted.

    Parameters
    ----------
    N : integer
        The size of the network (number of nodes).
    p : float between 0 and 1
        Probability of links.
    directed : Boolean
        True if a directed graph is desired, False if an undirected graph is
        desired.
    selfloops: Boolean
        True if self-loops are allowed, False otherwise.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN
        The adjacency matrix of the generated random graph.

    See Also
    --------
    RandomGraph : Random graphs with specified number of links
    """
    # 0) SECURITY CHECKS. Make sure that adequate parameters are given.
    if (p < 0.0 or p > 1.0):
        raise ValueError( "Probability p out of bounds. Insert value between 0 and 1" )

    # 1) FOR DIRECTED GRAPHS
    # 1.1) Create a random matrix with normally distributed values
    adjmatrix = numpy.random.rand(N,N)
    # 1.2 Select the entries with value <= p
    adjmatrix = np.where(adjmatrix <= p, 1, 0).astype(outdtype)

    # 2) FOR UNDIRECTED GRAPHS
    if not directed:
        adjmatrix[np.tril_indices(N)] = 0
        adjmatrix += adjmatrix.T

    # 3) Remove the diagonal if no self-loops are desired
    if not selfloops:
        adjmatrix[np.diag_indices(N)] = 0

    return adjmatrix

def RandomGraph(N, L, directed=False, selfloops=False):
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
                raise ValueError( "L out of bounds, max(L) = N**2 =", maxL )
        else:
            maxL = N*(N-1)
            if L > maxL:
                raise ValueError( "L out of bounds, max(L) = N*(N-1) =", maxL )
    else:
        if selfloops:
            maxL = 0.5*N*(N+1)
            if L > maxL:
                raise ValueError( "L out of bounds, max(L) = 1/2*N*(N+1) =", maxL )
        else:
            maxL = 0.5*N*(N-1)
            if L > maxL:
                raise ValueError( "L out of bounds, max(L) = 1/2*N*(N-1) =", maxL )

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

def BarabasiAlbertGraph(N, m, outdtype=np.uint8):
    """Generates scale-free networks after the Barabasi & Albert model.

    The Barabasi and Albert model (Science 286 (1999)) creates networks with
    a scale-free degree distribution of exponent gamma = -3 by a growth
    process with preferential attachment. At each iteration a new node is
    included that connects to the existing nodes with probability proportional
    to their degree.
    In our implementation the network is initialized by a complete graph of
    size m + 1 nodes.

    Parameters
    ----------
    N : integer
        Number of nodes of the final network.
    m : integer
        Number of links that each new node makes during the growing process.
    outdtype : numpy compatible data type, optional.
        The data type of the resulting adjacency matrix of the scale-free
        network.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN
        The adjacency matrix of the generated scale-free network.
    """

    # 1) INITIATE THE NETWORK AS A COMPLETE GRAPH OF SIZE m
    adjmatrix = np.zeros((N,N),outdtype)
    adjmatrix[:m+1,:m+1] = np.ones((m+1,m+1),outdtype)
    adjmatrix[np.diag_indices(m+1)] = 0

    # 2) PERFORM THE PREFERENTIAL ATTACHMENT GROWTH
    # 2.0) Create a list initially containing the hubs ki times
    nodelist = list(range(m+1))*m

    # 2.1) Add a new node
    for i in range(m+1,N):
        counter = 0
        neighbours = []
        while counter < m:
            # 2.2) Choose a node at random
            jid = int(len(nodelist)*numpy.random.rand())
            j = nodelist[jid]
            # Avoid double links
            if adjmatrix[i,j]: continue

            # 2.3) Connect the new node with the randomly chosen node
            adjmatrix[i,j] = 1
            adjmatrix[j,i] = 1
            neighbours.append(j)
            counter += 1

        # 2.4) Update nodelist for new probability in the next iteration
        nodelist += neighbours
        nodelist += m*[i]

    return adjmatrix

def ScaleFreeGraph(N, L, exponent=3.0, directed=False):
    """Generates scale-free graphs of given size and exponent.

    It follows the method proposed by Goh, Kahng & Kim 'Universal Behaviour
    of Load Distribution in SF networks' PRL 87 (2001). Every node is chosen
    with probability alpha = 1 / (gamma - 1) where gamma is the desired
    exponent.

    Parameters
    ----------
    N : integer
        Number of nodes of the final network.
    L : integer
        Number of links of the resulting random network.
    exponent : float (optional)
        The exponent (in positive) of the degree-distribution of the resulting
        networks. Recommended values between 2 and 3.
    directed : boolean, optional
        False if a graph is desired, True if digraphs are desired. In case
        of digraphs, both the input and the output degrees follow a scale-
        free distribution but uncorrelated between them.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN.
        The adjacency matrix of the generated scale-free network.

    Notes
    -----
    In case of directed networks the input and the output degrees of the
    nodes are correlated, e.g., input hubs are also output hubs.
    """
    # 0) SECURITY CHECKS
    if directed:
        maxL = N*(N-1)
        if L > maxL:
            raise ValueError( "L out of bounds, max(L) = N*(N-1) =", maxL )
    else:
        maxL = 0.5*N*(N-1)
        if L > maxL:
            raise ValueError( "L out of bounds, max(L) = 1/2*N*(N-1) =", maxL )

    # 1) PREPARE FOR THE CALCULATIONS
    adjmatrix = np.zeros((N,N), np.uint8)

    # Create a degree sequence
    alpha = 1.0/(exponent - 1.0)
    nodeweights = (np.arange(N) +1)**-alpha
    # Probability of a node to be chosen
    nodeweights /= nodeweights.sum()
    nodecumprobability = nodeweights.cumsum()
    del nodeweights

    # 2) CONNECT THE NETWORK
    counter = 1
    while counter <= L:

        # 2.1) Choose two nodes to connect
        xhead = numpy.random.rand()     # A random number between 0 and 1
        xsum = np.sum(np.sign(nodecumprobability-xhead))
        head = int(0.5*(N-xsum))

        xtail = numpy.random.rand()
        xsum = np.sum(np.sign(nodecumprobability-xtail))
        tail = int(0.5*(N-xsum))

        # 2.2) Do not allow self loops and multiple edges
        if head == tail: continue
        if adjmatrix[head,tail]: continue

        # 2.3) If conditions are satisfied, make the link
        adjmatrix[head,tail] = 1
        if not directed:
            adjmatrix[tail,head] = 1
        counter += 1

    return adjmatrix


############################################################################
"""NETWORK REWIRING ALGORITHMS"""
def RewireNetwork(adjmatrix, prewire=10, directed=None, weighted=False):
    """Randomises an input graph conserving the degrees of its nodes.

    It uses the link switching method to rewire networks while conserving
    both the input and the output degrees of the nodes. See references:
    A.R. Rao & S. Bandyopadhyay, Sankhya, Ser. A 58, 225 (1996), and
    J.M. Roberts, Soc. Networks 22, 273 (2000).

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network to be rewired. 'adjmatrix' itself
        won't be rewired but a new matrix is returned.
    prewire : float, optional
        Fraction of links to be rewired. See Usage for further instructions.
    directed : Boolean, optional
        Specifies the directedness of the returned network. If 'directed' is
        None (default), the function checks the directedness of the input
        network and performs the rewiring accordingly. To save computational
        time the parameter can be manually specified.
        If 'directed' is True, a directed network is returned regardless
        'adjmatrix' is directed or undirected. Input and output degrees are
        conserved.
        If 'directed' is False, an undirected network is returned. Only valid
        when 'adjmatrix' is undirected. DO NOT set 'directed = False' when
        'adjmatrix' is a directed network, degress won't be conserved and the
        function won't raise an error.
    weighted : Boolean, optional
        Specifies whether the weights of the links are conserved.
        If 'weighted' is False, a binary network is returned.
        If 'weighted' is True, links are switched conserving their weights.
        Notice, however, that it is impossible to conserve simultaneously
        both the degrees and the strength of all nodes. Hence, the weights of
        the links are conserved but only the input strength of the nodes
        are conserved. The output strength of nodes are not conserved.

    Returns
    -------
    rewmatrix : ndarray of rank-2
        The dtype of 'rewmatrix' depends on the option 'weighted'. If
        'weighted' is True, 'rewmatrix' is of same dtype as 'adjmatrix'.
        If 'weighted' is False, 'rewmatrix' has dtype uint8.

    Usage
    -----
    The parameter 'prewire' is not a probability, but a parameter that
    controls the number of iterations the link swithching procedure happens.
    At each iteration two links are switched, therefore for a given value of
    'prewire', prewire*L links are rewired in 1/2*prewire*L iterations, where
    L in the number of links. To make sure that all structure of the network
    has been completely randomized (appart from the degrees), values of
    prewire >= 5 are recommended.
    """
    # 0) PREPARE FOR THE CALCULATIONS
    # 0.1) Check the conditions for the rewiring process
    if directed==None:
        recip = Reciprocity(adjmatrix)
        if recip == 1.0: directed = False
        else: directed = True
    if weighted:
        rewmatrix = adjmatrix.copy()
    else:
        rewmatrix = np.where(adjmatrix,1,0).astype(np.uint8)

    N = len(rewmatrix)
    # 0.2) Generate the list of links
    if directed:
        linklist = np.array(rewmatrix.nonzero())
    else:
        # Apply nonzero only to the upper triangular part of the matrix
        linklist = np.array(np.triu(rewmatrix).nonzero())

    L = len(linklist[0])
    iterations = int(round(0.5*prewire*L))

    # DO THE REWIRING
    count = 0
    while count < iterations:
        # 1) SELECT TWO LINKS AT RANDOM:
        linkid1 = int(L*numpy.random.rand())
        linkid2 = int(L*numpy.random.rand())
        # Security check. If the two links are the same, discard the iteration
        if linkid1 == linkid2: continue

        h1 = linklist[0,linkid1]; t1 = linklist[1,linkid1]
        h2 = linklist[0,linkid2]; t2 = linklist[1,linkid2]

        # 2) SECURITY CHECKS TO AVOID INTRODUCING UNDESIRED LINKS
        # Avoid formation of self-loops
        if h1 == t2: continue
        if h2 == t1: continue
        # Avoid formation of double-links
        if rewmatrix[h1,t2]: continue
        if rewmatrix[h2,t1]: continue
        # Avoid trivial exchange of links
        if h1 == h2: continue
        if t1 == t2: continue

        # 3) IF ALL CONDITIONS SUCCESFUL, REWIRE
        # 3.1) Rewire the matrix
        if directed:
            # Put the new links
            if weighted:
                rewmatrix[h1,t2] = rewmatrix[h2,t2]
                rewmatrix[h2,t1] = rewmatrix[h1,t1]
            else:
                rewmatrix[h1,t2] = 1
                rewmatrix[h2,t1] = 1

            # Remove the old links
            rewmatrix[h1,t1] = 0
            rewmatrix[h2,t2] = 0

        else:
            # Put the new links
            if weighted:
                rewmatrix[h1,t2] = rewmatrix[h2,t2]
                rewmatrix[t2,h1] = rewmatrix[t1,h1]
                rewmatrix[h2,t1] = rewmatrix[h1,t1]
                rewmatrix[t1,h2] = rewmatrix[t2,h2]
            else:
                rewmatrix[h1,t2] = 1; rewmatrix[t2,h1] = 1
                rewmatrix[h2,t1] = 1; rewmatrix[t1,h2] = 1

            # Remove the old links
            rewmatrix[h1,t1] = 0.0; rewmatrix[t1,h1] = 0.0
            rewmatrix[h2,t2] = 0.0; rewmatrix[t2,h2] = 0.0

        # 3.2) Update the linklist
        linklist[1,linkid1] = t2
        linklist[1,linkid2] = t1

        # 3.3) Count the succesful realization
        count += 1

    return rewmatrix

def ModularityPreservingGraph(adjmatrix, partition, directed=None, selfloops=None):
    """Randomises an input graph conserving its modular structure.

    Given the adjacency matrix of a graph and a partition of its nodes, this
    function returns a network that contains the same number of links in
    each community and across them, but with the links randomly seeded.
    modules and across them as in the input graph.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network to be rewired. 'adjmatrix' itself
        won't be rewired but a new matrix is returned.
    partition : list of ndarrays of dtype = uint
        A list containing the indices of the nodes in each module.
    directed : Boolean (optional)
        True if a directed graph is desired, False if an undirected graph is
        desired.
    selfloops: Boolean (optional)
        True if self-loops are allowed, False otherwise.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and dtype = uin8
        The adjacency matrix of the generated random modular graph.

    See Also
    --------
    HMRandomGraph : Generates Nested random hierarchical/modular networks.
    RewireNetwork : Randomises an input graph conserving degrees of nodes.
    ModularHeterogenousGraph : Generates random modular networks of given
                               module sizes and densities.
    """

    # 0) SECURITY CHECKS AND SETUP
    # Convert the input network into a boolean matrix
    adjmatrix = adjmatrix.astype('bool')

    # Check if the original network is directed or undirected
    if directed == None:
        if Reciprocity(adjmatrix) == 1: directed = False
        else: directed = True

    # Check if the original network accepts self-loops
    if selfloops == None:
        if adjmatrix.trace(): selfloops = True
        else: selfloops = False

    # 1) INITIATE THE NEW ADJACENCY MATRIX AND HELPERS
    N = len(adjmatrix)
    ncoms = len(partition)
    randmatrix = np.zeros((N,N), np.uint8)

    # 2) GENERATE THE BLOCK-WISE RANDOM GRAPH
    for c1, com1 in enumerate(partition):
        N1 = len(com1)

        # 2.1) Seed the random links within the module
        submatrix = ExtractSubmatrix(adjmatrix, com1)
        Lblock = submatrix.astype('bool').sum()
        if not directed:
            Lblock = int( round(np.float32(Lblock / 2)) )
        counter = 0
        while counter < Lblock:
            # Pick up two nodes at random
            idx1 = int( N1*numpy.random.rand() )
            idx2 = int( N1*numpy.random.rand() )
            # Check for self-loops
            if not selfloops and idx1==idx2: continue
            source = com1[idx1]
            target = com1[idx2]

            # Check if they are already connected
            if randmatrix[source,target] == 1: continue

            # If the nodes are linkable, place the link
            randmatrix[source,target] = 1
            if not directed:
                randmatrix[target,source] = 1

            counter += 1

        # 2.2) Seed the random links across modules
        for c2 in range(ncoms):
            if not directed and c2 < c1: continue
            if c1 == c2: continue

            com2 = partition[c2]
            N2 = len(com2)
            subnet = ExtractSubmatrix(adjmatrix,com1,com2)
            Lblock = subnet.astype('bool').sum()
            counter = 0
            while counter < Lblock:
                # Pick up two nodes at random, each from a different module
                idx1 = int( N1*numpy.random.rand() )
                idx2 = int( N2*numpy.random.rand() )
                source = com1[idx1]
                target = com2[idx2]

                # Check if they can be linked, otherwise look for another pair
                if randmatrix[source,target] == 1: continue

                # If the nodes are linkable, place the link
                randmatrix[source,target] = 1
                if not directed:
                    randmatrix[target,source] = 1

                counter += 1

    return randmatrix


############################################################################
"""MODULAR AND HIERARCHICAL NETWORK MODELS"""
def ModularHeterogeneousGraph(Nsizelist, pintlist, pext, directed=False, selfloops=False):
    """
    Generates random modular networks of given module sizes and densities.

    This function generates modular networks in which both the internal
    links within the modules and the external links across modules are shed
    at random. As in the Erdos-Renyi model, every pair of nodes is linked
    with the specified probability and therefore, the total number of links
    slightly varies from one realization to another. The user will specify
    the size and the internal connection probability for each of the modules,
    and one unique probability for the external links across communities.

    Parameters
    ----------
    Nsizelist : list, tuple or array of integers
        A list containing the desired size (number of nodes) for every
        module in the network.
    pintlist : list, tuple or array of floats
        A list  containing the internal link probability for each of the
        modules. All values must range between 0 and 1.
    pext : float
        The external probability of connection between nodes in different
        communities.
    directed : Boolean (optional)
        True if a directed graph is desired, False if an undirected graph is
        desired.
    selfloops: Boolean (optional)
        True if self-loops are allowed, False otherwise.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and dtype = uin8
        The adjacency matrix of the generated random modular graph.
    partition : list of ndarrays of dtype = uint
        A list containing the indices of the nodes in each module.

    Usage and examples
    ------------------
    Setting Nsizelist = [100,200,300], pintlist = [0.3, 0.3, 0.5] and
    pext = 0.01 will generate a random graph of size N = 600 nodes with
    three modules of sizes N1 = 100, N2 = 200 and N3 = 300 nodes respectively.
    Each module is an Erdos-Renyi random graph with link probability
    (approximate density of links) pint1 = 0.3, pint2 = 0.3 and pint3 = 0.5
    respectively. The probability that a pair of nodes in different
    modules are connected is pext = 0.01.

    See Also
    --------
    HMRandomGraph : Generates Nested random hierarchical/modular networks.
    ErdosRenyiGraph : Generates random graphs with given link probability.
    """
    # 0) SECURITY CHECKS
    if len(Nsizelist) != len(pintlist):
        raise TypeError( "Parameters 'Nsizelist' and 'pintlist' not aligned." )
    if (pext < 0.0 or pext > 1.0):
        raise ValueError("Probability 'pext' out of bounds, insert value between 0 and 1" )
    for c in range(len(pintlist)):
        if (pintlist[c] < 0.0 or pintlist[c] > 1.0):
            raise ValueError( "Probability 'pintlist' out of bounds. Insert values between 0 and 1" )

    # 1) PREPARE TO CREATE THE NETWORK
    N = np.add.reduce(Nsizelist)
    ncommunities = len(Nsizelist)
    adjmatrix = np.zeros((N,N), np.uint8)

    # Define the partition
    counter = 0
    partition = []
    for c in range(ncommunities):
        com = np.arange(counter,counter+Nsizelist[c])
        partition.append(com)
        counter += Nsizelist[c]

    # 2) GENERATE THE RANDOM MODULAR NETWORK
    for c1, com1 in enumerate(partition):
        N1 = Nsizelist[c1]

        for c2, com2 in enumerate(partition):
            N2 = Nsizelist[c2]

            # 2.0) Choose the probability to use
            if c1 == c2: pthres = pintlist[c1]
            else: pthres = pext

            # 2.1) Create a random matrix with normally distributed values
            submatrix = numpy.random.rand(N1,N2)
            # 2.2) Select the entries with value <= p
            submatrix = np.where(submatrix <= pthres, 1, 0).astype(np.uint8)
            # 1.3) Copy the submatrix into the final adjmatrix
            imin = com1[0]; imax = com1[-1] + 1
            jmin = com2[0]; jmax = com2[-1] + 1
            adjmatrix[imin:imax,jmin:jmax] = submatrix

    # 3) IF GRAPH SHOULD BE UNDIRECTED...
    if not directed:
        adjmatrix[np.tril_indices(N, k=-1)] = 0
        adjmatrix += adjmatrix.T

    # 4) Remove the diagonal if no self-loops are desired
    if not selfloops:
        adjmatrix[np.diag_indices(N)] = 0

    return adjmatrix, partition

def HMpartition(HMshape):
    """Returns a partition of nodes for a hierarchical and modular network.
    network.

    Parameters
    ----------
    HMshape : list, tuple or array-like.
        HMshape must be a list of integer numbers indicating the target
        hierarchical and modular structure. For example, HMshape = [2,4,20]
        will create a network of 2 modules, each module is divided into 4
        modules of 20 nodes each.

    Returns
    -------
    partitions : list
        Contains the partition matrices for each of the nlevel-1 hierarchical
        levels of a network created with given parameter 'HMshape'.

    See Also
    --------
    HMRNetwork : Generates random modular and hierarchical networks.
    """
    N = np.multiply.reduce(HMshape)
    nlevels = len(HMshape)

    partitions = []
    for level in range(nlevels):
        # 2.1.1) Find the number of blocks per hierarchical level
        if level == 0: continue

        nblocks = int(np.multiply.reduce(HMshape[:level]))
        partitionmatrix = np.zeros((N,nblocks), np.uint8)

        # 2.2) For each block in the hierarchical level
        for b in range(nblocks):
            # 2.2.1) Find the indices for the block
            Nblock = N // nblocks
            iminblock = b*Nblock
            imaxblock = (b+1)*Nblock

            partitionmatrix[iminblock:imaxblock,b] = 1

        partitions.append(partitionmatrix)

    return partitions

def HMRandomGraph(HMshape, avklist, directed=False, outdtype=np.uint8):
    """Generates random hierarchical and modular networks of desired number
    of hierarchical levels and modules.

    This function generalizes the benchmark hierarchical and modular network
    models introduced in [M.E.J. Newman & M. Girvan, Phys. Rev. E 69, 026113
    (2004)] and in [A. Arenas, A. Diaz-Guilera & C.J. Perez-Vicente, Phys.
    Rev. Lett. 96, 114102 (2006)]. It creates networks of networks that are
    randomly connected at each hierarchical level, i.e., random networks
    that are randomly connected with each other forming larger modules that
    randomly connect with each other. The communities of each level are
    of same size and link density.

    Parameters
    ----------
    HMshape : list, tuple or array-like.
        HMshape must be a list of integer numbers indicating the target
        hierarchical and modular structure. For example, HMshape = [2,4,20]
        will create a network of 2 modules, each module is divided into 4
        modules of 20 nodes each.
    avklist : list, tuple or array-like.
        The mean degree of the nodes at each hierarchical level. For example,
        avklist = [1,3,20] will generate a network in which nodes will have,
        on average, 20 links with their neighbours at the lowest level
        submodule they belong to, 3 links with the rest of neighbours at the
        submodule they belong to at the second level, and 1 link with any
        node at any community of the rest of the network.
    directed : boolean, optional.
        If true, the resulting network will be directed, else, it will be
        undirected.
    outdtype : data type, optional
        Data-type of the resulting adjacency matrix. Default: uint8.

    Returns
    -------
    adjmatrix : ndarray of rank-2.
        The adjacency matrix of the generated modular and hierarchical
        network.

    Usage and examples
    ------------------
    The function accepts any desired number of hierarchical levels and
    partitions of the levels into desired number of nodes or communities.
    - HMshape = [4,4,16] will create a network of 4x4x16 = 256 nodes and
    three hierarchical levels. The network contains 4 communities, each
    subdivided into 4 communities of 16 nodes.
    - HMshape = [2,5,3,20] will generate a network of 2x5x3x20 = 600 nodes
    and four hierarchical levels. The network is composed of 2 communities,
    each subdivided into 5 communities, each subdivided into 3 communities
    of 20 nodes each.

    The parameter 'avklist' controls the density of connections between
    the modules in each of the hierarchical levels.
    - avklist = [1,3,13], for the first example above, means that each node
    is connected with 13 nodes of the first level, to 3 nodes in other
    communities at the second level and makes one link with nodes in the
    rest of the network which do not belong to the same modules in either
    the first or the second level.

    In order to generate the 13-4 or the 15-2 networks in [Arenas et al,
    Phys. Rev. Lett., 114102 (2006)], the function has to be called with
    parameters HMshape = [4,4,16] and avklist = [1,4,13] or avklist = [1,2,15]
    respectively.

    Notes
    -----
    The degrees given with parameter 'avklist' are average degrees. Due to
    the random generation process, every node will only have the desired
    internal and external degrees on the ensemble average, but not
    necessarily in each realization.

    See Also
    --------
    HMCentralisedGraph
        Random hierarchical and modular network with inter-modular
        connectivity centralised through local hubs.
    ModularHeterogeneousGraph
        Generates random modular networks with desired module sizes and
        densities.
    """
    def SeedLinks(adjmatrix, L, partition, directed=False):
        """This is a helper subfunction that seeds links at random.
        """
        # 0) SSECURITY CHECK
        if len(partition) < 2.0:
            raise TypeError( "Partition needs to have at least 2 communities" )

        ncoms = len(partition)
        Ncoms = len(partition[0])
        counter = 0
        while counter < L:
            # 1) Choose two communities at random
            com1idx = int(ncoms*numpy.random.rand())
            com2idx = int(ncoms*numpy.random.rand())
            # Make sure the nodes belong to different communities
            if com1idx == com2idx: continue

            # 2) Choose one node from each community
            nodeidx1 = int( Ncoms*numpy.random.rand() )
            nodeidx2 = int( Ncoms*numpy.random.rand() )
            node1 = partition[com1idx][nodeidx1]
            node2 = partition[com2idx][nodeidx2]
            # Avoid multiple links
            if adjmatrix[node1,node2]: continue

            # Link the nodes and update counter
            if directed:
                adjmatrix[node1,node2] = 1
                counter += 1
            else:
                adjmatrix[node1,node2] = 1
                adjmatrix[node2,node1] = 1
                counter += 1

    #______________________________________________________________________
    # 0) SECURITY CHECKS
    if len(HMshape) != len(avklist):
        raise ValueError( "HMshape and plist not aligned." )
    if HMshape[0] < 1:
        raise ValueError("HMshape[0] <= 1. First hierarchical level must contain more than one module." )

    # 1) PREPARE TO CREATE THE NETWORK
    N = np.multiply.reduce(HMshape)
    nlevels = len(HMshape)
    adjmatrix = np.zeros((N,N), outdtype)

    # 2) CREATE THE HM NETWORK BY SEEDING LINKS AT DIFFERENT SCALES
    # 2.1) For each hierarchical level
    for level in range(nlevels):
        # 2.1.1) Find the number of blocks per hierarchical level
        if HMshape[level] == 1: continue

        if level == 0: nblocks = 1
        else: nblocks = int(np.multiply.reduce(HMshape[:level]))

        # 2.2) For each block in the hierarchical level
        for b in range(nblocks):
            # 2.2.1) Find the indices for the block
            Nblock = N // nblocks
            iminblock = b*Nblock
            imaxblock = (b+1)*Nblock

            # 2.2.2) Create the partition of nodes into communities within
            # the current block in the current hierarchical level
            ncoms = HMshape[level]
            Ns = Nblock // ncoms
            partition = []
            for p in range(ncoms):
                partition.append(range(iminblock + p*Ns,iminblock + (p+1)*Ns))

            # 3) SEED THE LINKS AT RANDOM BETWEEN COMMUNITIES
            # NOTE: In the last hierarchical level, every node is one community
            if directed:
                Ls = Nblock * avklist[level]
            else:
                Ls = 0.5 * Nblock * avklist[level]

            SeedLinks(adjmatrix, Ls, partition, directed)

    return adjmatrix

def HMCentralisedGraph(HMshape, avklist, gammalist, directed=False, outdtype=np.uint8):
    """    Generates random hierarchical and modular networks of desired number
    of hierarchical levels and modules, with centralised inter-modular
    connectivity through local hubs.

    This function creates networks of networks that are randomly connected
    at each hierarchical level. The submodules at the lowest level are
    random scale-free-like graphs. Setting a high exponent for this
    level (e.g. gamma=[x,x,100]) will make them usual random graphs.
    At each heirarchical level, the links between modules are seeded
    following a preferential attachment rule (same as as for the generation
    of scale-free graphs). The nodes at every module are assigned different
    probability to link with other modules, hence, the inter-modular links
    end up concentrated on few hubs, with every community having its own
    set of hubs.

    For further details see reference:
    - G. Zamora-Lopez, Y. Chen et al. "Functional complexity emerging from
    anatomical constraints in the brain: the significance of network
    modularity and rich-clubs." Scientific Reports 6:38424 (2016).

    Parameters
    ----------
    HMshape : list, tuple or array-like.
        HMshape must be a list of integer numbers indicating the target
        hierarchical and modular structure. For example, HMshape = [4,50]
        will create a network composed of 4 modules of 50 nodes each.
        HMshape = [2,4,50] will create a network of 2 modules, each module
        divided into 4 submodules of 50 nodes each.
    avklist : list, tuple or array-like.
        The mean degree of the nodes at each hierarchical level. For example,
        avklist = [1,3,20] will generate a network in which nodes will have,
        on average, 20 links with their neighbours at the lowest level
        submodule they belong to, 3 links with the rest of neighbours at the
        submodule they belong to at the second level, and 1 link with any
        node at any community of the rest of the network.
    gammalist : list, tuple or array-like.
        Exponent controling the preferential attachement rule at each
        level. For a network with HMshape = [4,50], setting gammalist =
        [2.0,3.0] will create a network of four modules of 50 nodes each.
        Internally, the four submodules are scale-free-like networks with
        exponent gamma = 3.0. The links between modules are seeded at random
        but
    directed : boolean, optional.
        If true, the resulting network will be directed, else, it will be
        undirected.
    outdtype : data type, optional
        Data-type of the resulting adjacency matrix. Default: uint8.

    Returns
    -------
    adjmatrix : ndarray of rank-2.
        The adjacency matrix of the generated modular and hierarchical
        network.

    Usage and examples
    ------------------
    The function accepts any desired number of hierarchical levels and
    partitions of the levels into desired number of nodes or communities.
    - HMshape = [4,4,16] will create a network of 4x4x16 = 256 nodes and
    three hierarchical levels. The network contains 4 communities, each
    subdivided into 4 communities of 16 nodes.
    - HMshape = [2,5,3,20] will generate a network of 2x5x3x20 = 600 nodes
    and four hierarchical levels. The network is composed of 2 communities,
    each subdivided into 5 communities, each subdivided into 3 communities
    of 20 nodes each.

    The parameter 'avklist' controls the density of connections between
    the modules in each of the hierarchical levels.
    - avklist = [1,3,13], for the first example above, means that each node
    is connected with 13 nodes of the first level, to 3 nodes in other
    communities at the second level and makes one link with nodes in the
    rest of the network which do not belong to the same modules in either
    the first or the second level.

    The parameter 'gammalist' controls for the probability that every node
    has to be selected while links are seeded. This probability may change
    at different hierarchical levels.
    - gammalist = [2,100] with HMshape = [4,50], will lead to a modular
    network consisting of four random graphs (random because of the high
    gamma=100 at the level of the modules). To seed the inter-modular
    links, two nodes are chosen at random from two different communities.
    Internally, the probability of a node to be chosen differs such that
    each community will have its own hubs such that intermodular links
    preferentially occur between those hubs. The probability of a node within
    a module to be chosen as the target of an intermodular link is
    determined by the exponent gamma = 2, which would eventually to a
    schale-free graph of gamma = 2 in the limit of large networks.

    See Also
    --------
    HMRandomGraph
        Generates random hierarchical/modular network of given shape.
    ModularHeterogeneousGraph
        Generates random modular networks with desired module sizes and
        densities.
    """
    def GenerateBlock(L, partition, cumprobability, directed=False):
        """This is a helper function which generates the modules at a given
        hierachical level.
        """
        if len(partition) < 1:
            raise ValueError( "No modules given, Generating a random graph" )

        ncoms = len(partition)
        Ns = len(partition[0])
        N = ncoms*Ns

        blockmatrix = np.zeros((N,N),np.uint8)
        counter = 0
        while counter < L:
            # 1) Choose two communities at random
            com1idx = int(ncoms*numpy.random.rand())
            com2idx = int(ncoms*numpy.random.rand())
            # Make sure the communities are different
            if com1idx == com2idx: continue

            # 2) Choose two nodes to connect given the cumulative probabilities
            x = numpy.random.rand()     # A random number between 0 and 1
            xsum = np.sum(np.sign(cumprobability-x))
            idx = int(0.5*(Ns-xsum))
            head = com1idx * Ns + idx

            x = numpy.random.rand()
            xsum = np.sum(np.sign(cumprobability-x))
            idx = int(0.5*(Ns-xsum))
            tail = com2idx * Ns + idx

            # 3.2) Do not allow self loops and multiple edges
            if head == tail: continue
            if blockmatrix[head,tail]: continue

            # 3.3) If conditions are satisfied, make the link
            blockmatrix[head,tail] = 1
            if not directed:
                blockmatrix[tail,head] = 1
            counter += 1

        return blockmatrix

    def SkewedRandomGraph(L, cumprobability, directed=False):
        """This is a helper function to seed inter-modular links.
        """
        Ncom = len(cumprobability)

        blockmatrix = np.zeros((Ncom,Ncom), np.uint8)
        counter = 0
        while counter < L:
            # 1) Choose two nodes to connect given the cumulative probabilities
            x = numpy.random.rand()
            xsum = np.sum(np.sign(cumprobability-x))
            head = int(0.5*(Ncom-xsum))

            x = numpy.random.rand()
            xsum = np.sum(np.sign(cumprobability-x))
            tail = int(0.5*(Ncom-xsum))

            # 3.2) Do not allow self loops and multiple edges
            if head == tail: continue
            if blockmatrix[head,tail]: continue

            # 3.3) If conditions are satisfied, make the link
            blockmatrix[head,tail] = 1
            if not directed:
                blockmatrix[tail,head] = 1
            counter += 1

        return blockmatrix

    #______________________________________________________________________
    # 0) SECURITY CHECKS
    if len(HMshape) != len(avklist):
        raise ValueError( "HMshape and avklist are not alighned." )
    if gammalist:
        if len(gammalist) != len(HMshape):
            raise ValueError( "HMshape and gammalist are not alighned." )

    # 1) PREPARE TO CREATE THE NETWORK
    N = np.multiply.reduce(HMshape)
    nlevels = len(HMshape)
    adjmatrix = np.zeros((N,N), outdtype)

    # 1.1) If no hub parameters given, connect modules at random.
    # This case returns the same networks as 'HMRandomGraph' function.
    if not gammalist:
        alpha = np.ones(nlevels, float)
    else:
        alpha = 1.0 / (np.array(gammalist,float) - 1.0)

    # 2) CREATE THE HM NETWORK BY SEEDING LINKS AT DIFFERENT SCALES
    for level in range(nlevels-1):
        # 2.1) Find the number of blocks, communities and nodes
        # Number of blocks
        if HMshape[level] == 1: continue
        if level == 0: nblocks = 1
        else:
            nblocks = int(np.multiply.reduce(HMshape[:level]))
        # Number of nodes per block
        Nblock = np.multiply.reduce(HMshape[level:])
        # Number of communities per block
        ncoms = HMshape[level]
        # Number of nodes per community
        Ncom = Nblock // ncoms

        # 2.2) Define typical partition for the current hierarchical level
        partition = np.zeros((ncoms,Ncom), np.uint)
        for i in range(ncoms):
            partition[i] = np.arange(i*Ncom,(i+1)*Ncom, dtype=np.uint)

        # 2.3) Define the typical node selection probabilities within a block
        if level < nlevels - 2:
            nodeweights = np.ones(Ncom,np.float)
            ncomsnext = HMshape[level+1]
            Ncomnext = Ncom // ncomsnext
            for i in range(ncomsnext):
                nodeweights[i*Ncomnext:(i*Ncomnext+Ncomnext)] = \
                    ((np.arange(Ncomnext)+1).astype(float))**-alpha[level]
            # Probability of a node to be chosen
            nodeweights /= nodeweights.sum()
            cumprobability = nodeweights.cumsum()
        else:
            nodeweights = np.ones(Ncom,np.float)
            nodeweights = ((np.arange(Ncom)+1).astype(float))**-alpha[level]
            # Probability of a node to be chosen
            nodeweights /= nodeweights.sum()
            cumprobability = nodeweights.cumsum()

        # 2.4) Determine the number of links to be seed per block
        if directed:
            Lblock = Nblock * avklist[level]
        else:
            Lblock = 0.5 * Nblock * avklist[level]

        # 2.5) Seed intercommunity links in the current hierarchical level
        for b in range(nblocks):
            minidx = b*Nblock
            maxidx = (b+1)*Nblock
            adjmatrix[minidx:maxidx,minidx:maxidx] = \
                GenerateBlock(Lblock, partition, cumprobability, directed)

    # 3) CREATE THE LAST HIERARCHICAL LEVEL
    # 3.1 ) Find the number communities and nodes per community
    Ncom = HMshape[-1]
    ncoms = len(adjmatrix) // Ncom

    # 3.2) Define the typical node selection probabilities within a community
    nodeweights = ((np.arange(Ncom)+1).astype(np.float))**(-alpha[-1])
    nodeweights /= nodeweights.sum()    # Probability of a node to be chosen
    cumprobability = nodeweights.cumsum()

    # 3.3) Determine the number of links to be seed per community
    if directed:
        Lcom = Ncom* avklist[-1]
    else:
        Lcom = 0.5 * Ncom * avklist[-1]

    # 3.4) Seed intercommunity links in the current hierarchical level
    for i in range(ncoms):
        minidx = i*Ncom
        maxidx = (i+1)*Ncom

        adjmatrix[minidx:maxidx,minidx:maxidx] = \
            SkewedRandomGraph(Lcom, cumprobability, directed)

    return adjmatrix

def RavaszBarabasiGraph(Nmotif=4, hlevels=3, hublinks=True):
    """Generates a Ravasz-Barabasi hierarchical graph.

    The Ravasz-Barabasi graph is a pseudo-fractal network that starts from a
    basic motif, a complete graph of size 'Nmotif'. Further hierarchical levels
    are constructed replicating the original motif, by making each of the nodes
    in original motif as the central node of the new replicated motif. All new
    'peripheral' vertices link

    The network is claimed to produce scale-free-like degree distribution but
    modular structure is rather unclear. Two different versions of the model
    where published. In Ref. [science2002], the hubs of each hierarchical level
    are connected forming a ring graph. Note, however, that these links
    between the hubs do not lead to a rich-club. The model in Ref. [pre2003],
    ignored the connections among the local hubs at same hierarchical level hubs.
    This function creates both versions, controled by parameter 'hublinks'.

    Parameters
    ----------
    Nmotif : integer
        Size of the basic motif.
    hlevels : integer
        Number of hierarchical levels. Notice that hlevels = 1 will return the
        basic motif, a complete graph of size Nmotif.

    Returns
    -------
    adjmatrix : ndarray of rank-2.
        The adjacency matrix of the generated scale-free network. The size of
        the final network is N = Nmotif^hlevels.

    References
    ----------
    [science2002] E. Ravasz, A.L. Somera, D.A. Mongru, Z.N. Oltvai and A.L.
    Barabasi "Hierarchical organization of modularity in metabolic networks."
    Science 297, 1551-1555 (2002).
    [pre2003] E. Ravasz and A. Barabasi "Hierarchical organization in complex
    networks" Phys. Rev. E 67, 026112 (2003).
    """
    # 0) SECURITY CHECKS
    if Nmotif < 4:
        raise ValueError('Basic motif needs at least four nodes, Nmotif > 3.')
    if hlevels < 1:
        raise ValueError('Number of hirarchical levels needs to be positive, hlevels > 0.')

    # The number of nodes of one module in three different scales
    N = Nmotif**(hlevels)
    adjmatrix = np.zeros((N,N), np.uint8)

    # 1) DEFINE THE BASIC MOTIF
    adjmatrix[:Nmotif,:Nmotif] = 1
    adjmatrix[np.diag_indices(Nmotif)] = 0

    # 2) GROW AND CONNECT THE NETWORK
    for lev in range(1,hlevels):
        hublist = np.arange(Nmotif**lev,Nmotif**(lev+1),Nmotif**lev)

        # 2.1) Link the hubs at the new level forming a ring
        if hublinks:
            for i in range(Nmotif-2):
                node1 = hublist[i]
                node2 = hublist[i+1]
                adjmatrix[node1,node2] = 1
                adjmatrix[node2,node1] = 1
            adjmatrix[hublist[0],hublist[-1]] = 1
            adjmatrix[hublist[-1],hublist[0]] = 1

        # 2.2) Define the basis motif for the current level
        motif = adjmatrix[:Nmotif**lev,:Nmotif**lev]
        newNmotif = len(motif)

        # 2.3) Make copies of the motif at each hub location
        for hub in hublist:
            adjmatrix[hub:hub+newNmotif,hub:hub+newNmotif] = motif

        # 2.4) Connect all new peripheral nodes with the central hub
        for i in range(Nmotif**lev,Nmotif**(lev+1)):
            if i not in hublist:
                adjmatrix[0,i] = 1
                adjmatrix[i,0] = 1

    return adjmatrix

#
