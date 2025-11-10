# -*- coding: utf-8 -*-
# Copyright (c) 2018, Gorka Zamora-López and Romain Brasselet
# <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
=============================================
ULTRA-SHORT AND ULTRA-LONG NETWORK GENERATORS
=============================================

This module contains functions to generate ultra-short and ultra-long networks,
that is, networks with shortest possible pathlength (smallest efficiency) and
longest possible pathlength (largest efficiency) of arbitrary number of nodes N
and number of links L. Generators for both graphs and digraphs are given.

For ultra-long networks, in some cases no apparent model exist to generate the
networks for each individual L. In those cases, solutions for special cases are
provided. These are the solutions of order M < N, which add L_M = 1/2 M(M-1)
links simultaneously.

ULTRA-SHORT NETWORKS
====================
USgraph
    Generates an ultra-short graph of size N and number of edges L.
USgraph_Disconnected
    Generates a disconnected graph with largest possible efficiency.
USgraph_Random
    Generates a connected ultra-short graph, random case.
USgraph_RichClub
    Generates a connected ultra-short graph, rich-club configuration.

USdigraph
    Generates an ultra-short directed graph of size N and number of edges L.
USdigraph_iDirectedRing
    Generates an incomplete directed ring.
USdigraph_iStarDigraph
    Generates an incomplete star digraph
USdigraph_FlowerDigraph
    Generates a flower digraph.
USdigraph_Random
    Generates a connected ultra-short digraph, random configuration.
USdigraph_RichClub
    Generates a connected ultra-short digraph, rich-club configuration.

ULTRA-LONG NETWORKS
====================
ULgraph_Connected
    Generates a connected ultra-long graph of specified N and L.
ULgraph_Disconnected_Mcomplete
    Generates an M-complete graph of specified N and order M.

ULdigraph_Connected_Range1_MBS
    Generates a connected ultra-long directed graph of order M.
ULdigraph_Connected_Range2
    Generates a connected ultra-long digraph of density rho >= 1/2 + 1/N.
ULdigraph_Disconnected_Range1
    Generates a digraph with smallest possible efficiency, L <= 1/2 N(N-1).
ULdigraph_Disconnected_Range2
    Generates a digraph with smallest possible efficiency, when L >= 1/2 N(N-1)
    and order M.

"""
# Standard library imports
import warnings
# Third party imports
import numpy as np
import numpy.random
# Local imports
from .limits import Effic_iDirectedRing, Effic_iStarDigraph, Effic_FlowerDigraph


## ULTRA-SHORT NETWORKS ########################################################
## UNDIRECTED GRAPHS __________________________________________________________
def USgraph(N,L, uscase='Random', reportdisco=False):
    """
    Generates an ultra-short graph of size N and number of edges L.

    - If L < N-1, returns a disconnected star graph, with infinite pathlength
    but largest efficiency possible.
    - If L = N-1, returns a start graph.
    - If L > N-1, the ultra-short theorem states that any graph with diameter
    diam(G) = 2 is ultra-short, regardless of the precise organization of its
    links. Among all possible ultra-short graphs, the function returns graphs
    following two optional models (parameter 'uscase').
    usecase='Random' returns a star graph with additional edges seeded randomly.
    usecase='RichClub' returns a star graph with additional edges are seeded
    such that consecutive hubs are created. Leading to ultra-short graphs with
    a rich-club.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.
    uscase : string (optional)
        A string with the value 'Random' or 'RichClub'. It specifies the type
        of ultra-short graph to generate when L > 2(N-1).
    reportdisco : boolean (optional)
        If 'True', informs that the generated graph is disconnected.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    USgraph_Disconnected :
    USgraph_Random :
    USgraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2:  raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:  raise ValueError( "L does not take negative values." )
    Ltot = 0.5*N*(N-1)
    if L > int(Ltot):  raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )
    assert uscase=='Random' or uscase=='RichClub', "'uscase' not valid. Please enter either 'Random' or 'RichClub'."

    # 1) SELECT THE MODEL AND GENERATE THE NETWORK
    if L < N-1:
        # Generate a disconnected US graph
        adjmatrix = USgraph_Disconnected(N,L)
        if reportdisco:
            warnings.warn("Generated digraph is disconnected", category=RuntimeWarning)
    else:
        # Generate a connected US graph
        if uscase == 'Random':
            adjmatrix = USgraph_Random(N,L)
        elif uscase == 'RichClub':
            adjmatrix = USgraph_RichClub(N,L)

    return adjmatrix

def USgraph_Disconnected(N,L):
    """
    Generates a disconnected graph with largest possible efficiency.

    When L < N-1, graphs are necessarily disconnected (thus have infinite
    pathlength). This function creates a disconnected star graph, i.e, a graph
    with largest possible efficiency in the range 0 < L < N-1.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    USgraph :
    USgraph_Random :
    USgraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:    raise ValueError( "L does not take negative values" )
    if L > N-1:  raise ValueError( "L out of bounds. max(L) = (N-1)" )

    # 1) CREATE THE INCOMPLETE STAR GRAPH
    adjmatrix = np.zeros((N,N), np.uint8)
    adjmatrix[0,1:L+1] = 1
    adjmatrix[1:L+1,0] = 1

    return adjmatrix

def USgraph_Random(N,L):
    """
    Generates a connected ultra-short graph, random case.

    Generates a graph with shortest possible pathlength (largest possible
    efficiency) by adding edges randomly to an initial star graph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    USgraph :
    USgraph_Disconnected :
    USgraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < N-1:       raise ValueError( "L out of bounds, min(L) = N-1" )
    if L > int(Ltot): raise ValueError( "L out of bounds, max(L) = 1/2 * N*(N-1)" )

    # 1) CREATE THE INITIAL STAR GRAPH
    adjmatrix = np.zeros((N,N), np.uint8)
    adjmatrix[0,1:] = 1
    adjmatrix[1:,0] = 1

    # 2) SEED THE REMAINING ARCS AT RANDOM
    counter = N-1
    while counter < L:
        # 2.1) Pick up two nodes at random
        source = int(N * numpy.random.rand())
        target = int(N * numpy.random.rand())

        # 2.2) Check if they can be linked, otherwise look for another pair
        if adjmatrix[source,target] == 1: continue
        if source == target: continue

        # 2.3) If the nodes are linkable, place the link
        adjmatrix[source,target] = 1
        adjmatrix[target,source] = 1
        counter += 1

    return adjmatrix

def USgraph_RichClub(N,L):
    """
    Generates a connected ultra-short graph, rich-club case.

    Generates a graph with shortest possible pathlength (largest possible
    efficiency) by seeding links to the initial star graph orderly, which leads
    to a rich-club architecture, e.g., an overlap of several star graphs with
    the central hubs interconnected.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    USgraph :
    USgraph_Disconnected :
    USgraph_Random :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < N-1:       raise ValueError( "L out of bounds, min(L) = N - 1" )
    if L > int(Ltot): raise ValueError( "L out of bounds, max(L) = 1/2 * N*(N-1)" )

    # 1) CREATE THE INITIAL STAR GRAPH
    adjmatrix = np.zeros((N,N), np.uint8)
    adjmatrix[0,1:] = 1
    adjmatrix[1:,0] = 1

    # 2) SEED THE REMAINING ARCS, FORMING A RICH-CLUB
    counter = N-1
    # 2.1) Stop if the network is a star graph ( L = 2(N-1) )
    if counter == L:
        finished = True
    elif counter < L:
        finished = False

    # 2.2) Else, seed the remaining edges favouring the formation of hubs.
    for i in range(1,N):
        if finished: break
        for j in range(i+1,N):
            # Place the edge
            adjmatrix[i,j] = 1
            adjmatrix[j,i] = 1
            # Check if job is finished
            counter += 1
            if counter == L:
                finished = True
                break

    return adjmatrix


## DIRECTED GRAPHS ____________________________________________________________
def USdigraph(N,L, onlyconnected=True, uscase='Random', reportdisco=False):
    """
    Generates an ultra-short directed graph of size N and number of edges L.

    Ultra-short digraph generation involves up to five different digraph models.
    This function sorts the models based on two optional parameters
    'onlyconnected' and 'uscase'.

    If L < N, any digraph is necessarily disconnected. Here, two models compete
    for the largest efficiency: Incomplete Directed Rings and Incomplete Star
    Digraphs. The function will identify which of the two has the largest
    efficiency, based on N and L, and return a digraph of the winner model.

    If N <= L < 2(N-1), optimal digraphs may be either connected or disconnected.
    If only connected digraphs are desired (onlyconnected=True), then it will
    return a Flower Digraph. But if digraphs are allowed to be disconnected
    (onlyconnected=False), then Flower Digraphs compete with Incomplete Star
    Digraphs for the largest efficiency. The function will identify which of
    the two has the largest efficiency, based on N and L, and return a digraph
    of the winner model.

    If L > 2(N-1) arcs, the ultra-short theorem states that any directed graph
    with diam(G) = 2 is ultra-short, regardless of its precise organization.
    Among all possible ultra-short graphs, the function returns digraphs
    following two optional models (parameter 'uscase').
    usecase='Random' returns a star graph with additional arcs seeded randomly.
    usecase='RichClub' returns a star graph with additional arcs seeded such
    that consecutive hubs are created. Leading to ultra-short graphs with
    a rich-club.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    L : integer
        Number of directed arcs of the digraph.
    onlyconnected : boolean (optional)
        If true, function will return Flower Digraphs when N < L < 2(N-1).
        If false, then the function will evaluate whether a Flower Digraph or an
        incomplete Directed Ring has the largest efficiency and return the
        winner model.
    uscase : string (optional)
        A string with the value 'Random' or 'RichClub'. It specifies the type
        of ultra-short digraph to generate in the dense case, when L > 2(N-1).
    reportdisco : boolean (optional)
        If 'True', basic feedback will be reported about the generated digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    USdigraph_iDirectedRing :
    USdigraph_iStarDigraph :
    USdigraph_FlowerDigraph :
    USdigraph_Random :
    USdigraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:   raise ValueError( "L does not take negative values." )
    Ltot = N*(N-1)
    if L > int(Ltot):  raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )
    assert uscase=='Random' or uscase=='RichClub', "'uscase' not valid. Please enter either 'Random' or 'RichClub'."

    # 1) CHOOSE CORRECT MODEL AND GENERATE THE CORRESPONDING ULTRA-SHORT DIGRAPH
    if L < N:
        # Digraph is disconnected. Identify the optimal model.
        # Disconnected directed rings and disconnected Star Digraphs compete
        effic1 = Effic_iDirectedRing(N,L)
        effic2 = Effic_iStarDigraph(N,L)

        # Generate the digraph with shortest efficiency
        if effic1 >= effic2:
            adjmatrix = USdigraph_iDirectedRing(N,L)
        else:
            adjmatrix = USdigraph_iStarDigraph(N,L)
        if reportdisco:
            warnings.warn("Generated digraph is disconnected", category=RuntimeWarning)

    elif L < 2*(N-1):
        # Digraph may be disconnected.
        if onlyconnected:
            # Only connected digraphs wanted
            adjmatrix = USdigraph_FlowerDigraph(N,L)
        else:
            # Disconnected digraphs allowed (we want largest efficiency)
            # Flower digraphs and disconnected Star Digraphs compete
            effic1 = Effic_FlowerDigraph(N,L)
            effic2 = Effic_iStarDigraph(N,L)

            # Generate the digraph with largest efficiency
            if effic1 >= effic2:
                adjmatrix = USdigraph_FlowerDigraph(N,L)
            else:
                adjmatrix = USdigraph_iStarDigraph(N,L)
                if reportdisco:
                    warnings.warn("Generated digraph is disconnected", category=RuntimeWarning)

    else:
        # The digraph is connected
        if uscase == 'Random':
            adjmatrix = USdigraph_Random(N,L)
        elif uscase == 'RichClub':
            adjmatrix = USdigraph_RichClub(N,L)

    return adjmatrix

def USdigraph_iDirectedRing(N,L):
    """
    Generates an incomplete directed ring.

    This model is restricted to the cases when L < N. An incomplete directed
    ring consists of a directed ring of size N' = L, with the remaining N - L
    nodes being isolated.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    L : integer
        Number of directed arcs of the digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    USdigraph :
    USdigraph_iStarDigraph :
    USdigraph_FlowerDigraph :
    USdigraph_Random :
    USdigraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:  raise ValueError( "L does not take negative values." )
    if L > N:  raise ValueError( "L out of bounds. max(L) = N" )

    # 1) CREATE THE NETWORK, CONTAINING ONE CYCLE OF LENGTH L
    adjmatrix = np.zeros((N,N), np.uint8)

    if L < 2:
        nodes = range(L+1)
    else:
        nodes = [0] + list(range(1,L)) + [0]

    # Seed the arcs and create the cycle
    for i in range(L):
        node1 = nodes[i]
        node2 = nodes[i+1]
        adjmatrix[node1,node2] = 1

    return adjmatrix

def USdigraph_iStarDigraph(N,L):
    """
    Generates an incomplete star digraph

    This model is restricted to L <= 2(N-1). An incomplete star digraph consists
    of a star graph of size N' = L / 2 <= N with the remaining N - N' nodes
    being isolated. If L is even, then the hub is bidirectionally connected to
    L/2 nodes. If L is odd, then the hub is bidirectionally connected to
    floor(L/2) nodes and the remaining arc points from the hub to one of the
    isolated nodes.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    L : integer
        Number of directed arcs of the digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    USdigraph :
    USdigraph_iDirectedRing :
    USdigraph_FlowerDigraph :
    USdigraph_Random :
    USdigraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:        raise ValueError( "L does not take negative values." )
    if L > 2*(N-1):  raise ValueError( "L out of bounds. max(L) = 2*(N-1)" )

    # 1) CREATE THE NETWORK, A STAR DIGRAPH
    adjmatrix = np.zeros((N,N), np.uint8)

    # 1.1) Seed the links
    Lund = int(L/2)
    adjmatrix[0,1:Lund+1] = 1
    adjmatrix[1:Lund+1,0] = 1

    # 1.2) Add an extra link if L is an odd number
    if 2*Lund < L:
        adjmatrix[0,Lund+1] = 1

    return adjmatrix

def USdigraph_FlowerDigraph(N,L):
    """
    Generates a flower digraph.

    This model is restricted to the range N <= L <= 2(N-1). A flower digraph
    consists of a set of directed rings (the petals) which overlap in a single
    hub. If possible, All petals should be of same size, otherwise, the smallest
    petal is only one node sorter than the longest. Special cases of flower
    digraphs:
    When L = N, the flower digraph is a directed ring.
    When L = 2(N-1), the flower digraph is a star graph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    L : integer
        Number of directed arcs of the digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    USdigraph :
    USdigraph_iDirectedRing :
    USdigraph_iStarDigraph :
    USdigraph_Random :
    USdigraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < N:     raise ValueError( "L out of bounds, min(L) = N." )
    if L > 2*(N-1): raise ValueError( "L out of bounds. max(L) = 2*(N-1)" )

    # 1) CALCULATE THE NUMBER OF CYCLES AND THEIR SIZES
    nloops = L - (N-1)
    minloopsize = int(L / nloops)
    Nlooplist = minloopsize * np.ones(nloops, np.uint64)

    # Make some cycles longer (+1) if needed
    difference = L - int(Nlooplist.sum())
    if difference > 0:
        for i in range(difference):
            Nlooplist[i] += 1
        # Sort by sizes
        Nlooplist = Nlooplist[::-1]

    # 2) CREATE THE NETWORK BY ADDING THE CYCLES
    adjmatrix = np.zeros((N,N), np.uint8)

    counter = 1
    for c in range(nloops):
        # Create the cycle as a list of nodes
        Nloop = int(Nlooplist[c])
        nodes = [0] + list(range(counter, counter+Nloop-1)) + [0]
        counter += (Nloop-1)

        # Seed the arcs in the current cycle:
        for i in range(Nloop):
            node1 = nodes[i]
            node2 = nodes[i+1]
            adjmatrix[node1,node2] = 1

    return adjmatrix

def USdigraph_Random(N,L):
    """
    Generates a connected ultra-short digraph, random configuration.

    Generates a digraph with shortest possible pathlength (largest possible
    efficiency) by adding arcs randomly to an initial star graph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the idgraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    USdigraph :
    USdigraph_iDirectedRing :
    USdigraph_iStarDigraph :
    USdigraph_FlowerDigraph :
    USdigraph_RichClub :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lmin = 2*(N-1); Ltot = N*(N-1)
    if L < Lmin:  raise ValueError( "L out of bounds, min(L) = 2*(N-1)" )
    if L > Ltot:  raise ValueError( "L out of bounds. max(L) = N*(N-1)" )

    # 1) CREATE THE INITIAL STAR GRAPH
    adjmatrix = np.zeros((N,N), np.uint8)
    adjmatrix[0,1:] = 1
    adjmatrix[1:,0] = 1

    # 2) SEED THE REMAINING ARCS AT RANDOM
    counter = 2*(N-1)
    while counter < L:
        # 2.1) Pick up two nodes at random
        source = int(N * numpy.random.rand())
        target = int(N * numpy.random.rand())

        # 2.2) Check if they can be linked, otherwise look for another pair
        if adjmatrix[source,target] == 1: continue
        if source == target: continue

        # 2.3) If the nodes are linkable, place the link
        adjmatrix[source,target] = 1
        counter += 1

    return adjmatrix

def USdigraph_RichClub(N,L):
    """
    Generates a connected ultra-short digraph, rich-club configuration.

    Generates a digraph with shortest possible pathlength (largest possible
    efficiency) by seeding links to an initial star graph orderly, which leads
    to a rich-club architecture, i.e., the overlap of several star graphs with
    the central hubs interconnected.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of directed arcs in the graph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    USdigraph :
    USdigraph_iDirectedRing :
    USdigraph_iStarDigraph :
    USdigraph_FlowerDigraph :
    USdigraph_Random :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lmin = 2*(N-1); Ltot = N*(N-1)
    if L < Lmin:  raise ValueError( "L out of bounds, min(L) = 2*(N-1)" )
    if L > Ltot:  raise ValueError( "L out of bounds. max(L) = N*(N-1)" )

    # 1) CREATE THE INITIAL STAR GRAPH
    adjmatrix = np.zeros((N,N), np.uint8)
    adjmatrix[0,1:] = 1
    adjmatrix[1:,0] = 1

    # 2) SEED THE REMAINING ARCS, FORMING A RICH-CLUB
    counter = 2*(N-1)

    # 2.1) Stop if the network is a star graph ( L = 2(N-1) )
    if counter == L:
        finished = True
    elif counter < L:
        finished = False

    # 2.2) Else, seed the remaining arcs favouring the formation of hubs.
    for i in range(1,N):
        if finished: break
        for j in range(i+1,N):
            # Place the arc from hub to leaf, and check if job is finished
            adjmatrix[i,j] = 1
            counter += 1
            if counter == L:
                finished = True
                break
            # Place the reciprocal arc (leaf to hub) and check if job is finished
            adjmatrix[j,i] = 1
            counter += 1
            if counter == L:
                finished = True
                break

    return adjmatrix


## ULTRA-LONG NETWORKS ########################################################
## UNDIRECTED GRAPHS __________________________________________________________
def ULgraph_Connected(N,L):
    """
    Generates a connected ultra-long graph of specified N and L.

    This function returns CONNECTED ultra-long graphs, which consist of an
    orderly addition of edges to an initial path graph (or line graph).
    These graphs have the longest possible pathlength and the smallest
    efficiency, only if, connected graphs are considered. In the case of the
    ultra-long boundary, the networks with smallest efficiency are always
    disconnected but such configurations are not always unique. See
    function ULgraph_Mcomplete(N,M) for a generator of special cases, of graphs
    (disconnected) with smallest possible efficiency (E = density).

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)
    D. Barnpoutis & R.M. Murray *Extremal properties of complex networks*
    arXiv:1104.5532 (2011).
    L. Gulyas, G. Horvath, T. Cseri & G. Kampis *An estimation of the
    shortest and largest average path length in graphs of given density*
    arXiv:1101.2549 (2011).

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    ULgraph_Mcomplete :
    """
    # 0) SECURITY CHECK
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lmin = (N-1); Ltot = 0.5*N*(N-1)
    if L < Lmin:  raise ValueError( "L out of bounds, min(L) = (N-1)" )
    if L > Ltot:  raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )

    # 1) GENERATE THE NETWORK
    # Generate the initial path graph
    adjmatrix = np.eye(N,k=1,dtype=np.uint8)

    if L > N-1:
        # Seed the remaining edges orderly
        counter = N-1
        finished = False
        for j in range(2,N):
            if finished: break
            for i in range(j-1):
                adjmatrix[i,j] = 1
                counter += 1
                if counter >= L:
                    finished = True
                    break

    adjmatrix += adjmatrix.T
    return adjmatrix

def ULgraph_Disconnected_Mcomplete(N,M):
    """
    Generates an M-complete graph of specified N and order M.

    In the case of the ultra-long boundary, the networks with smallest efficiency
    are always disconnected but such configurations are not always unique.
    This function return graphs with smallest efficiency (DISCONNECTED ultra-
    long digraphs) for specific values of L = 1/2 M(M-1), with M <= N. These
    consist of a complete graph of size M and N-M isolated nodes. The
    efficiency of these special cases is (E = density).

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    M : integer
            Size of the embedded complete subgraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated graph.

    See Also
    --------
    ULgraph_Connected :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if M < 0:  raise ValueError( "M does not take negative values" )
    if M > N:  raise ValueError( "M out of range, max(M) = N" )

    # 1) GENERATE THE NETWORK
    # Initiate the empty network
    adjmatrix = np.zeros((N,N), np.uint8)

    # Add the complete subgraph of order M
    triuidx = np.triu_indices(M,k=1)
    adjmatrix[triuidx] = 1
    adjmatrix += adjmatrix.T

    return adjmatrix


## DIRECTED GRAPHS ____________________________________________________________
def ULdigraph_Connected_Range1_MBS(N,M):
    """
    Generates a connected ultra-long directed graph of order M.

    The generation of strongly connected ultra-long digraphs undergoes a
    transition when L = (N-1) + 1/2 N(N-1), or density rho = 1/2 + 1/N.
    If L < (N-1) + 1/2 N(N-1) algorithmic solutions exists only for particular
    values of L. These digraphs consist of a directed ring with the remaining
    arcs seeded following the opposite orientation, i.e., each of the first M
    nodes point to all its predecessors. We name such arcset as an M-Backwards
    Subgraph, with M <= N.
    These digraphs have longest possible pathlength and smallest efficiency,
    if only if, strongly connected digraphs are considered. The digraphs with
    smallest efficiency are always disconnected.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    M : integer
        Order of the M-Backward Subgraph. M = 2, 3, 4, ... , N-1.
        If M = 1, returns a directed ring.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated ultra-long digraph.

    See Also
    --------
    ULdigraph_Connected_Intermediate :
    ULdigraph_Connected_Range2 :
    """

    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if M < 1:  raise ValueError( "M out of range, min(M) = 1" )
    if M > N-1:  raise ValueError( "M out of range, max(M) = N-1" )

    # 1) GENERATE THE NETWORK
    # Create the baseline directed ring
    adjmatrix = np.zeros((N,N), np.uint8)
    for i in range(N-1):
        adjmatrix[i,i+1] = 1
    adjmatrix[-1,0] = 1

    # Seed the arcs for the M-backward subgraph
    for i in range(M):
        for j in range(i):
            adjmatrix[i,j] = 1

    return adjmatrix

def ULdigraph_Connected_Intermediate(N):
    """
    Connected ultra-long directed graph with L = (N-1) + 1/2 N(N-1) arcs.

    The generation of strongly connected ultra-long digraphs undergoes a
    transition when L = (N-1) + 1/2 N(N-1), or density rho = 1/2 + 1/N.
    If L = (N-1) + 1/2 N(N-1) a unique configuration exists with longest
    pathlength. This consists of the superposition of a directed ring with a
    complete directed acyclic graph whose arcs are pointing in the opposite
    orientations as the ring. That is, each vertex j points to all other
    vertices with i < j.
    These digraphs have longest possible pathlength and smallest efficiency,
    if only if, strongly connected digraphs are considered. The digraphs with
    smallest efficiency are always disconnected.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated ultra-long digraph.

    See Also
    --------
    ULdigraph_Connected_Range1_MBS :
    ULdigraph_Connected_Range2 :
    """

    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )

    # 1) GENERATE THE NETWORK
    adjmatrix = np.ones((N,N), np.uint8)
    triuidx = np.triu_indices(N,k=2)
    adjmatrix[triuidx] = 0
    adjmatrix[np.diag_indices(N)] = 0

    return adjmatrix

def ULdigraph_Connected_Range2(N,L):
    """
    Generates a connected ultra-long digraph of density rho >= 1/2 + 1/N.

    The generation of strongly connected ultra-long digraphs undergoes a
    transition when L = (N-1) + 1/2 N(N-1), or density rho = 1/2 + 1/N.
    At this point, the ultra-long digraph corresponds to a directed ring with
    the largest possible M-backward added.
    If L > (N-1) + 1/2 N(N-1) further arcs are added in the same orientation as
    the directed ring, and bilateralising the arcs of the M-Backwards Subgraph,
    leading to complete subgraphs within the network.
    These digraphs have longest possible pathlength and smallest efficiency,
    if only if, strongly connected digraphs are considered. The digraphs with
    smallest efficiency are always disconnected.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of arcs in the digraph. L >= 1/2*N*(N-1) + (N-1)

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated ultra-long digraph.

    See Also
    --------
    ULdigraph_Connected_Range1_MBS :
    ULdigraph_Connected_Intermediate :
    """

    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lstart = int( 0.5*N*(N-1) + (N-1) );    Ltot = N*(N-1)
    if L < Lstart:  raise ValueError( "L out of range, min(L) = 1/2*N*(N-1) + (N-1)" )
    if L > Ltot:    raise ValueError( "L out of range, max(L) = N*(N-1)." )

    # 1) CREATE THE INITIAL NETWORK (DENSEST M-BACKWARDS SUBGRAPH)
    adjmatrix = ULdigraph_Connected_Intermediate(N)

    # 2) SEED THE REMAINING FORWARD ARCS ORDERLY
    Lr = L - Lstart
    if Lr:
        counter = 0
        finished = False
        for j in range(2,N):
            if finished: break
            for i in range(j-1):
                if adjmatrix[i,j]: continue
                adjmatrix[i,j] = 1
                counter += 1
                if counter >= Lr:
                    finished = True
                    break

    return adjmatrix

def ULdigraph_Disconnected_Range1(N,L):
    """
    Generates a digraph with smallest possible efficiency, L <= 1/2 N(N-1).

    Disconected digraphs with largest possible efficiency are represented by
    two different models. Both lead to disconnected digraphs.
    If L < 1/2 N(N-1) then the optimal configuration is a directed acyclic
    graph consisting of adding an M-Backwards subgraphs to an initially empty
    network. In this case arcs can be seeded one-by-one.
    If L = N(N-1), the procedure creates the densest possible DAG, which we
    refer to as the Complete DAG.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph. L <= 1/2 N(N-1)

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    ULdigraph_Disconnected_Range2 :
    """
    # 0) SECURITY CHECK
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lmax = int( 0.5*N*(N-1) )
    if L < 0:    raise ValueError( "L does not take negative values" )
    if L > Lmax: raise ValueError( "L out of range, max(L) = 0.5 * N*(N-1)." )

    # 1) CREATE THE NETWORK
    adjmatrix = np.zeros((N,N), np.uint8)

    # Seed the edges
    # NOTE. Fix this! I don't like to launch the algorithm this way, "manually"
    # for L > 0. The algorithm should be written in a manner that for L = 0,
    # it naturally stops withouth adding any arc.
    if L > 0:
        counter = 0
        finished = False
        for i in range(1,N):
            if finished: break
            for j in range(i):
                adjmatrix[i,j] = 1
                counter += 1
                if counter == L:
                    finished = True
                    break

    return adjmatrix

def ULdigraph_Disconnected_Range2(N,M):
    """
    Generates a digraph with smallest possible efficiency, when L >= 1/2 N(N-1)
    and order M.

    Disconected digraphs with largest possible efficiency are represented by
    two different models. Both lead to disconnected digraphs.
    If L = N(N-1), the optimal configuration is the densest possible DAG, which
    we refer to as the Complete DAG.
    If L > 1/2 N(N-1), digraphs with smalles efficiency are generated by orderly
    bilateralising the arcs of a Complete DAG. Unique solutions exist for
    given values of adding groups of L_M = 1/2 M(M-1) arcs (with M <=N) in the
    forward direction such that the first M nodes receive arcs from all their
    predecessors. We name such arcsets as M-Forward Subgraphs. These digraphs
    all have efficiency E = density.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    M : integer
        Order of the M-Forward Subgraph. M = 2, 3, 4, ... , N.

    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the generated digraph.

    See Also
    --------
    ULdigraph_Diconnected_Range1 :
    """
    # 0) SECURITY CHECK
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if M < 2:  raise ValueError( "M out of range, min(M) = 2" )
    if M > N:  raise ValueError( "M out of range. max(M) = N" )

    # 1) GENERATE THE NETWORK
    # Create the initial complete DAG
    adjmatrix = np.zeros((N,N), np.uint8)
    trilidx = np.tril_indices(N,k=-1)
    adjmatrix[trilidx] = 1

    # Seed the arcs for the M-forward subgraph
    for j in range(1,M):
        for i in range(j):
            adjmatrix[i,j] = 1

    return adjmatrix



##
