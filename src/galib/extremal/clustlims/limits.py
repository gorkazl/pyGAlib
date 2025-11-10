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
================================================
ANALYTIC BOUNDARIES OF PATHLENGTH AND EFFICIENCY
================================================

This module contains functions to calculate, analytically, the upper and the
lower limits for the pathlength and efficiency of networks with arbitrary
number of nodes N and number of links L. Results are available for
graphs and digraphs, in separate functions.

ULTRA-SHORT LIMITS
==================
Pathlen_USgraph
    Shortest possible average pathlength for graphs (ultra-short limit).
Effic_USgraph
    Largest possible global efficiency for graphs (ultra-short limit).

Pathlen_USdigraph
    Shortest possible average pathlength for directed graphs (ultra-short limit).
Pathlen_FlowerDigraph
    Average pathlength of a flower digraph (ultra-short limit).
Effic_USdigraph
    Largest possible global efficiency for digraphs (ultra-short limit).
Effic_FlowerDigraph
    Global efficiency of a flower digraph (ultra-short limit).
Effic_iDirectedRing
    Global efficiency of an incomplete directed ring (ultra-short limit).
Effic_iStarDigraph
    Global efficiency of an incomplete star digraph (ultra-short limit).

ULTRA-LONG LIMITS
=================
Pathlen_ULgraph
    Longest possible average pathlength for graphs (ultra-long limit).
Effic_ULgraph
    Smallest possible efficiency for graphs (ultra-long limit).

Pathlen_ULdigraph
    Longest possible pathlength for strongly connected digraphs (ultra-long limit).
Pathlen_ULdigraph_Range1_MBS
    Longest possible pathlength of a strongly connected digraphs containing an
    M-Backwards subgraph (ultra-long limit).
    Valid in the range when L <= (N-1) + 1/2 N(N-1), or density <= 1/2 + 1/N.
Pathlen_ULdigraph_Range1_Approx
    Longest possible pathlength for a strongly connected digraph (ultra-long
    limit, approximation).
    Valid in the range when L <= (N-1) + 1/2 N(N-1), or density <= 1/2 + 1/N.
Pathlen_ULdigraph_Intermediate
    Largest possible pathlength for a digraph with L = (N-1) + 1/2 N(N-1) arcs.
Pathlen_ULdigraph_Range2
    Longest possible pathlength for strongly connected digraphs (ultra-long limit).
    Valid in the range when L >= (N-1) + 1/2 N(N-1), or density >= 1/2 + 1/N.
Effic_ULdigraph_Range1_MBS:
    Smallest possible efficiency for strongly connected digraphs (ultra-long limit)
    containing an M-Backwards subgraph.
    Valid in the range when L <= (N-1) + 1/2 N(N-1), or density <= 1/2 + 1/N.
Effic_ULdigraph_Intermediate
    Smallest possible efficiency for strongly connected digraphs (ultra-long limit)
    containing exactly L = (N-1) + 1/2 N(N-1) arcs, or density = 1/2 + 1/N.
Effic_ULdigraph_Disconnected
    Smallest possible efficiency of a digraph that is not strongly connected.

"""
# Standard library imports
import warnings
# Third party imports
import numpy as np
import numpy.random
import scipy.special
# Local imports

# Global constants
EMconstant = 0.57721566490153286060


## ULTRA-SHORT LIMIT ##########################################################
## UNDIRECTED GRAPHS __________________________________________________________
def Pathlen_USgraph(N,L):
    """
    Shortest possible average pathlength for graphs (ultra-short limit).

    Calculates the shortest possible pathlength a graph of arbitrary number of
    nodes N and number of undirected edges L could possibly have. In case that
    L < N-1, the graph is necessarily disconnected and the function returns an
    infinite value.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)
    For the connected case, see also:
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
    avpathlen : float
        The average pathlength of an undirected ultra-short graph.

    See Also
    --------
    Effic_USgraph :
    USgraph :
    USdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )

    # 1) CALCULATE THE PATHLENGTH
    if L < N-1:
        # Network is disconnected, infinite pathlength
        avpathlen = np.inf
    else:
        # Network is dense enough for at least one connected configuration
        avpathlen = 2.0 - L / Ltot

    return avpathlen

def Effic_USgraph(N,L):
    """
    Largest possible global efficiency for graphs (ultra-short limit).

    Calculates the largest possible efficiency that a graph of arbitrary number
    of nodes N and number of undirected edges L could possibly have. Networks
    with shortest possible pathlength have largest efficiency. In case
    L < N-1, the graph is necessarily disconnected and the function returns the
    efficiency of a disconnected star graph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)
    For the connected case, see also
    D. Barnpoutis & R.M. Murray *Extremal properties of complex networks*
    arXiv:1104.5532 (2011).

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.

    Returns
    -------
    efficiency : float
        The efficiency of an undirected ultra-short graph.

    See Also
    --------
    Pathlen_USgraph :
    USgraph :
    USdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )

    # 1) CHOOSE CORRECT MODEL AND CALCULATE CORRESPONDING EFFICIENCY
    if L < N:
        # Network is disconnected, return efficiency of incomplete star graph
        efficiency = 0.25 * (L**2 + 3*L) / Ltot
    else:
        # Network is connected and dense enough to form a star graph (diam = 2)
        efficiency = 0.5 * ( 1.0 + L / Ltot )

    return efficiency

## DIRECTED GRAPHS __________________________________________________________
def Pathlen_USdigraph(N,L):
    """
    Shortest possible average pathlength for directed graphs (ultra-short limit).

    Calculates the shortest possible pathlength a connected directed graph of
    arbitrary number of nodes N and number of directed arcs L could possibly
    have. In case N < L < 2(N-1), its returns the pathlength of flower digraphs.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of a strongly connected ultra-short digraph.

    See Also
    --------
    Pathlen_FlowerDigraph :
    Effic_USdigraph :
    USgraph :
    USdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds. max(L) = N*(N-1)" )

    # 1) CHOOSE CORRECT MODEL AND CALCULATE THE CORRESPONDING PATHLENGTH
    if L < N:
        # Network is disconnected, infinite pathlength
        avpathlen = np.inf
    elif L < 2*(N-1):
        # Network is sparse but connected. Flower digraphs.
        avpathlen = Pathlen_FlowerDigraph(N,L)
    else:
        # Network is connected and dense enough to form a star graph (diam = 2)
        avpathlen = 2.0 - float(L) / Ltot

    return avpathlen

def Pathlen_FlowerDigraph(N,L):
    """
    Average pathlength of a flower digraph (ultra-short limit).

    When the number of arcs is in the range N < L < 2(N-1), strongly connected
    digraphs exist but they are not dense enough for the ultra-short theorem to
    apply. Their diameter is dima(G) > 2. Flower digraphs are the special
    configurations which minimise the pathlength in this regime. This function
    thus returns the pathlength of a flower digraph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of a flower digraph.

    See Also
    --------
    Effic_USdigraph :
    USgraph :
    USdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < N:       raise ValueError( "L out of range, min(L) = N" )
    if L > 2*(N-1): raise ValueError( "L out of range, max(L) = 2*(N-1)" )

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

    # 2) CALCULATE THE PATHLENGTH ANALYTICALLY
    # 2.1) Find the number of loops of size x and x+1
    x = Nlooplist.min()
    mx = len(np.where(Nlooplist==x)[0])
    my = len(np.where(Nlooplist==x+1)[0])

    # 2.2) Calculate the contributions to pathlength per loop
    Ix = 0.5 * x**2 * (x-1)
    Iy = 0.5 * (x+1)**2 * x
    Jxx = np.float64( x*(x-1)**2 )
    Jyy = np.float64( x**2 * (x+1) )
    Jxy = np.float64( 0.5 * x*(x-1)*(2*x+1) )

    # 2.3) Calculate the terms of the final sum
    Dx = mx * (Ix + (mx-1) * Jxx)       # Loops of length d
    Dy = my * (Iy + (my-1) * Jyy)       # Loops of length d+1
    Dxy = 2.0 * mx*my * Jxy             # Between loops of length d and d+1

    Ltot = N*(N-1)
    avpathlen = ( Dx + Dy + Dxy ) / Ltot
    return avpathlen

def Effic_USdigraph(N,L):
    """
    Largest possible global efficiency for digraphs (ultra-short limit).

    Calculates the largest possible efficiency that a directed graph of
    arbitrary number of nodes N and number of directed arcs L could possibly
    have. Networks with shortest possible pathlength have largest efficiency.

    When L < N digraphs are disconnected and two models compete for largest
    efficiency, disconnected directed rings and disconnected star digraphs. The
    winning model depends on precise N and L.

    In the range N < L < 2(N-1), other two models compete: disconnected star
    digraph and flower digraphs. The winning model depends on precise N and L.
    In the range L >= 2(N-1), then the ultra-short theorem applies.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of a strongly connected ultra-short digraph.

    See Also
    --------
    Effic_FlowerDigraph :
    Effic_iDirectedRing :
    Effic_iStarDigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds. max(L) = N*(N-1)" )

    # 1) CHOOSE CORRECT MODEL AND CALCULATE THE CORRESPONDING PATHLENGTH
    if L <= N:
        # Network is disconnected, for sure
        # Disconnected directed ring and disconnected Star Digraphs compete
        effic1 = Effic_iDirectedRing(N,L)
        effic2 = Effic_iStarDigraph(N,L)
        efficiency = max(effic1, effic2)
    elif L < 2*(N-1):
        # Network is sparse. Could be connected or disconnected.
        # Flower digraphs and disconnected star digraphs compete
        effic1 = Effic_FlowerDigraph(N,L)
        effic2 = Effic_iStarDigraph(N,L)
        efficiency = max(effic1, effic2)
    else:
        # Network is connected and dense enough to form a star digraph
        # Ultra-short theorem applies
        efficiency = 0.5 * ( 1.0 + L / Ltot )

    return efficiency

def Effic_FlowerDigraph(N,L):
    """
    Global efficiency of a flower digraph (ultra-short limit).

    When the number of arcs is in the range N <= L <= 2(N-1), strongly connected
    digraphs exist but they are not dense enough for the ultra-short theorem to
    apply. Their diameter is diam(G) > 2. Flower digraphs are the special
    configurations which minimise the efficiency in this regime. This function
    thus returns the pathlength of a flower digraph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes of the digraph.
    L : integer
        Number of directed arcs of the digraph. For flower digraphs L must fall
        in the range N <= L <= 2(N-1).

    Returns
    -------
    efficiency : float
        The efficiency of the directed ultra-short graph.

    See Also
    --------
    Effic_USdigraph :
    Effic_iDirectedRing :
    Effic_iStarDigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2:   raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < N:   raise ValueError( "L out of range, min(L) = N" )
    if L > 2*(N-1): raise ValueError( "L out of range, max(L) = 2*(N-1)" )

    def Ix(x):
        result = x * (scipy.special.psi(x) + EMconstant)
        return result

    def Jxy(x,y):
        psi_xy = scipy.special.psi(x+y)
        psi_x = scipy.special.psi(x)
        psi_y = scipy.special.psi(y)
        result = (x+y-1) * psi_xy - x * psi_x - y * psi_y - (EMconstant+1)
        return result

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

    # 2) CALCULATE THE PATHLENGTH ANALYTICALLY
    # 2.1) Find the number of loops of size d and d+1
    x = Nlooplist.min()
    mx = len(np.where(Nlooplist==x)[0])
    my = len(np.where(Nlooplist==(x+1))[0])

    # 2.2) Calculate the loop contributions
    Dx = Ix(x)
    Dy = Ix(x+1)
    Dxx = Jxy(x,x)
    Dyy = Jxy(x+1,x+1)
    Dxy = Jxy(x,x+1)

    # 2.3) Calculate the terms of the final sum
    Ex = mx * (Dx + (mx-1) * Dxx)       # Loops of length x
    Ey = my * (Dy + (my-1) * Dyy)       # Loops of length x+1
    Exy = 2.0 * mx*my * Dxy             # Between loops of length d and x+1

    Ltot = N*(N-1)
    efficiency = (Ex + Ey + Exy) / Ltot

    return efficiency

def Effic_iDirectedRing(N,L):
    """
    Global efficiency of an incomplete directed ring (ultra-short limit).

    When the number of arcs is L < N, digraphs are necessarily disconnected.
    digraphs exist but they are not dense enough for the ultra-short theorem to
    apply. Disconnected directed rings are a special which maximise the
    efficiency in this regime. These consist of a directed ring of N' = L nodes
    and (N - N') isolated vertices. This function thus returns the efficiency of
    such a digraph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph. (L <= N must be given).

    Returns
    -------
    efficiency : float.
        The efficiency of an incomplete directed ring.

    See Also
    --------
    Effic_USdigraph :
    Effic_iStarDigraph :
    Effic_FlowerDigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:   raise ValueError( "L does not take negative values" )
    if L > N:   raise ValueError( "L out of range, max(L) = N" )

    # 1) CALCULATE THE EFFICIENCY
    Ltot = float(N*(N-1))
    if L == 0:
        efficiency = 0.0
    elif L == 1:
        efficiency = 1.0 / Ltot
    else:
        density = L / Ltot
        efficiency = density * (scipy.special.psi(L) + EMconstant)

    return efficiency

def Effic_iStarDigraph(N,L):
    """
    Global efficiency of an incomplete star digraph (ultra-short limit).

    When the number of arcs is L < 2(N-1), three models compete for largest
    efficiency. Incomplete star digraphs consist of a star graph of size N' < N,
    with remaining nodes isolated. In case L is 'odd' then one arc links the
    central hub with one of the isolated vertices. This function returns the
    efficiency of such a digraph.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph. L <= 2(N-1) must be given.

    Returns
    -------
    efficiency : float.
        The efficiency of an incomplete star digraph.

    See Also
    --------
    Effic_USdigraph :
    Effic_iDirectedRing :
    Effic_FlowerDigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2:  raise ValueError( "Network needs at least two nodes, N > 1" )
    if L < 0:  raise ValueError( "L does not take negative values" )
    if L > 2*(N-1): raise ValueError( "L out of range, max(L) = 2*(N-1)" )

    # 1) CALCULATE THE EFFICIENCY
    if L == 0:
        efficiency = 0.0
    else:
        Ltot = float(N*(N-1))
        undL = int(L/2)
        efficiency = (1.0/Ltot) * ( L + 0.5*undL * (L-undL-1) )

    return efficiency


## ULTRA-LONG LIMIT ###########################################################
## UNDIRECTED GRAPHS __________________________________________________________
def Pathlen_ULgraph(N,L):
    """
    Longest possible average pathlength for graphs (ultra-long limit).

    Calculates the longest possible pathlength a graph of arbitrary number of
    nodes N and number of undirected edges L could possibly have. In case that
    L < N-1, the graph is necessarily disconnected and the function returns an
    infinite value.

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
    avpathlen : float
        The average pathlength of an undirected ultra-long graph.

    See Also
    --------
    Effic_ULgraph :
    ULgraph :
    ULdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds. max(L) = 1/2 * N*(N-1)" )

    # 1) CALCULATE THE PATHLENGTH
    if L < N-1:
        # Network is disconnected, infinite pathlength
        avpathlen = np.inf
    else:
        # Number of nodes in the complete subgraph and in the tail
        Nc = 0.5 * (3 + np.sqrt(9 + 8*(L-N)))
        Nc = np.floor(Nc)
        Nt = N - Nc
        # Number of edges in complete subgraph
        Lc = 0.5 * Nc * (Nc-1)

        # Sum the contributions
        dijsum = Lc  - Nt * (L-N) + 1./6 * (N**3 - Nc**3 - 7*Nt)
        avpathlen = dijsum / Ltot

    return avpathlen

def Effic_ULgraph(N,L, connected=True):
    """
    Smallest possible efficiency for graphs (ultra-long limit).

    Calculates the largest possible efficiency that a graph of arbitrary number
    of nodes N and number of undirected edges L could possibly have. The outcome
    of optimal efficiency is splitted in two versions.

    If only the properties of connected graphs are desired, then it
    corresponds to the efficiency of the connected ultra-long graph, see also
    the function Pathlen_ULgraph(N,L).
    However, the smallest efficiency possible corresponds always to graphs that
    are disconnected. Given graphs are allowed to be disconnected, the smallest
    efficiency any graph could take equals the density of the graph.
    This choice is controlled by parameter 'connected' in the function.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)
    D. Barnpoutis & R.M. Murray *Extremal properties of complex networks*
    arXiv:1104.5532 (2011).

    Parameters
    ----------
    N : integer
        Number of nodes in the graph.
    L : integer
        Number of undirected edges in the graph.
    connected : boolean (optional)
        If 'True', the function returns the smallest efficiency given that the
        graph is connected. If 'False' it returns the absolute smallest
        efficiency, which corresponds to disconnected graphs.

    Returns
    -------
    efficiency : float
        The smallest efficiency of graph could take.

    See Also
    --------
    USgraph :
    ULgraph :
    Pathlen_USgraph :
    Effic_USgraph :
    Effic_ULgraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = 0.5*N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds, max(L) = 1/2 * N*(N-1)" )

    # 1) CHOOSE CORRECT MODEL AND CALCULATE CORRESPONDING EFFICIENCY
    if connected:
        # Allow only CONNECTED graphs
        if L < N-1:
            msg = "Number of edges out of range. L < N-1 incompatible with connected graphs. Returning NaN."
            warnings.warn(msg, RuntimeWarning)
            efficiency = np.nan
        else:
            # Number of nodes in the complete subgraph and in the tail
            Nc = 0.5 * (3 + np.sqrt(9 + 8*(L-N)))
            Nc = np.floor(Nc)
            Nt = N - Nc

            # Number of edges in complete subgraph and in the tail
            Lc = 0.5 * Nc * (Nc-1)
            Lt = Nt - 1
            Le = L - Lt - Lc

            # Sum the contributions
            eijsum = L - Nt - (Le-1)/(Nt +1) + \
                        N * (scipy.special.psi(Nt+2) + EMconstant - 1)
            efficiency = eijsum / Ltot

    else:
        # Allow efficiency of DISCONNECTED graphs
        efficiency = L / Ltot

    return efficiency


## DIRECTED GRAPHS ____________________________________________________________
def Pathlen_ULdigraph(N,L):
    """
    Longest possible pathlength for strongly connected digraphs (ultra-long limit).

    Calculates the longest possible pathlength a connected directed graph of
    arbitrary number of nodes N and number of directed arcs L could possibly
    have. The solutions are splitted at L' = 1/2 N(N-1) + (N-1).
    In the range N <= L <= L', exact solution is only known for special cases
    in which L = N + 1/2 M(M-1) with M < N. Here, we provide an approximation
    that can be estimated for any value of L.
    When L > L', then exact solutions are known.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of a strongly connected ultra-short digraph.

    See Also
    --------
    Pathlen_ULdigraph_Range1_MBS :
    Pathlen_ULdigraph_Range1_Approx :
    Pathlen_ULdigraph_Range2 :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of bounds, max(L) = N*(N-1)" )

    # 1) CHOOSE CORRECT MODEL AND CALCULATE THE CORRESPONDING PATHLENGTH
    Linterm = int( 0.5*N*(N-1) + (N-1) )

    if L < N:
        # Network is disconnected, infinite pathlength
        avpathlen = np.inf
    elif L < Linterm:
        # Number of arcs below the transition. Provide approximated result
        avpathlen = Pathlen_ULdigraph_Range1_Approx(N,L)
    elif L == Linterm:
        # Number of arcs at the transition.
        avpathlen = float(N+4) / 6
    else:
        # Number of arcs above the transition (densest regime).
        avpathlen = Pathlen_ULdigraph_Range2(N,L)

    return avpathlen

def Pathlen_ULdigraph_Range1_MBS(N,M):
    """
    Longest possible pathlength of a strongly connected digraphs containing an
    M-Backwards subgraph (ultra-long limit).
    Valid in the range when L <= (N-1) + 1/2 Lo, or density <= 1/2 + 1/N.

    Calculates the pathlength of a digraph consisting of a directed ring and
    and an M-backwards subgraph. An M-BS is a set of arcs where the first M
    nodes all point to all their predecessors. Such digraphs have
    L = N + 1/2 M(M-1) arcs and longest possible average pathlength, where
    M < N.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    M : integer
        Order of the M-Backward Subgraph. M = 2, 3, ... , N-1.
        If M = 1, returns the pathlength of a directed ring.

    Returns
    -------
    avpathlen : float
        The average pathlength of the digraph.

    See Also
    --------
    Pathlen_ULdigraph :
    Pathlen_ULdigraph_Range1_Approx :
    Pathlen_ULdigraph_Intermediate :
    Pathlen_ULdigraph_Range2 :
    """

    # 0) SECURITY CHECK
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if M < 1: raise ValueError( "M out of range, min(M) = 1" )
    if M > N-1: raise ValueError( "M out of range, max(M) = N-1." )

    # 1) CALCULATE THE PATHLENGTH
    M = float(M)
    term1 = 0.5*N
    term2 = 0.5*float(M*(M-1)) / (N*(N-1))
    term3 = N - float(M+4) / 3

    avpathlen = term1 - term2 * term3
    return avpathlen

def Pathlen_ULdigraph_Range1_Approx(N,L):
    """
    Longest possible pathlength for a strongly connected digraph (ultra-long
    limit, approximation).
    Valid in the range when L <= (N-1) + 1/2 N(N-1), or density <= 1/2 + 1/N.

    The solutions leading to digraphs with longest pathlength are splitted at
    L' = 1/2 N(N-1) + (N-1). In the range N <= L <= L', exact solution is only
    known for special cases in which L = N + 1/2 M(M-1) with M < N. This
    function returns an approximation of that pathlength that can be estimated
    for any value of L.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The approximated average pathlength of an ultra-short digraph.

    See Also
    --------
    Pathlen_ULdigraph :
    Pathlen_ULdigraph_Range1_MBS :
    Pathlen_ULdigraph_Intermediate :
    Pathlen_ULdigraph_Range2 :
    """
    # 0) Security checks
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lmin = N;   Lmax = 0.5*N*(N-1) + (N-1)
    if L < Lmin:      raise ValueError( "L out of range, min(L) = N" )
    if L > int(Lmax): raise ValueError( "L out of range, max(L) = 1/2*N*(N-1) + (N-1)" )

    # 1) Prepare for calculations
    Ltot = N*(N-1)
    density = float(L) / Ltot

    term1 = 1.0 + 3./2 * density
    term2 = (1.0 + 0.5*density) * np.sqrt(0.5*density)
    term3 = N * ( 0.5 - density + 1./3 * density * np.sqrt(2.0*density) )

    avpathlen = term1 - term2 + term3
    return avpathlen

def Pathlen_ULdigraph_Intermediate(N):
    """
    Largest possible pathlength for a digraph with L = (N-1) + 1/2 N(N-1) arcs.

    Calculates the pathlength of a digraph consisting of the superposition of
    a directed ring with a complete directed acyclic graph whose arcs are
    pointing in the opposite orientations as the ring. That is, each vertex j
    points to all other vertices with i < j.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of the digraph.

    See Also
    --------
    Pathlen_ULdigraph :
    Pathlen_ULdigraph_Range1_MBS :
    Pathlen_ULdigraph_Range1_Approx :
    Pathlen_ULdigraph_Range2 :
    """

    # 0) SECURITY CHECK
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )

    # 1) CALCULATE THE PATHLENGTH
    avpathlen = float(N + 4) / 6

    return avpathlen

def Pathlen_ULdigraph_Range2(N,L):
    """
    Longest possible pathlength of a strongly connected digraph (ultra-long limit).
    Valid in the range when L >= (N-1) + 1/2 N(N-1), or density >= 1/2 + 1/N.

    Calculates the longest possible pathlength a connected directed graph of
    arbitrary number of nodes N and number of directed arcs L could possibly
    have, whenever L > 1/2 N(N-1) + (N-1).

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of directed arcs in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of a strongly connected ultra-short digraph.

    See Also
    --------
    Pathlen_ULdigraph :
    Pathlen_ULdigraph_Range1_MBS :
    Pathlen_ULdigraph_Range1_Approx :
    Pathlen_ULdigraph_Intermediate :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Lstart = int( 0.5*N*(N-1) + (N-1) );    Ltot = N*(N-1)
    if L < Lstart:  raise ValueError( "L out of range, min(L) = 1/2*N*(N-1) + (N-1)" )
    if L > Ltot:    raise ValueError( "L out of range, max(L) = N*(N-1)." )

    # 1) CALCULATE THE CORRESPONDING PATHLENGTH
    Lr = L - Lstart
    term1 = float(N + 4.0)/6
    term2 = float(Lr) / N

    # Calculate term #3, which includes a series without a closed formula.
    partsum = 0
    finished = False
    counter = 0
    for i in range(1,N):
        if finished: break
        for j in range(i):
            if counter == Lr:
                finished = True
                break
            partsum += i
            counter += 1
    term3 = float(partsum) / Ltot

    # Sum all the contributiona
    avpathlen = term1 - term2 + term3
    return avpathlen

def Effic_ULdigraph_Range1_MBS(N,M):
    """
    Smallest possible efficiency for strongly connected digraphs (ultra-long limit)
    containing an M-Backwards subgraph.
    Valid in the range when L <= (N-1) + 1/2 N(N-1), or density <= 1/2 + 1/N.

    Calculates the efficiency of a digraph consisting of a directed ring and
    and an M-backwards subgraph. An M-BS is a set of arcs where the first M
    nodes all point to all their predecessors. Such digraphs have
    L = N + 1/2 M(M-1) arcs and longest possible average pathlength, where
    M < N.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    M : integer
        Order of the M-Backward Subgraph. M = 2, 3, ... , N-1.
        If M = 1, returns the pathlength of a directed ring.

    Returns
    -------
    avpathlen : float
        The average pathlength of the digraph.

    See Also
    --------
    Effic_ULdigraph_Intermediate :
    Effic_ULdigraph_Range2 :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    if M < 1: raise ValueError( "M out of range, min(M) = 1" )
    if M > N-1: raise ValueError( "M out of range, max(M) = N-1." )

    # 1) CALCULATE THE EFFICIENCY
    Lmax = N*(N-1)
    term1 = (scipy.special.psi(N) + EMconstant) / (N-1)
    term2 = float(M-1) / Lmax
    term3 = 0.5*M - scipy.special.psi(N)
    term4 = 0.0
    for j in range(1,M):
        term4 += scipy.special.psi(N-j)
    term4 = float(term4) / Lmax

    efficiency = term1 + term2 * term3 + term4
    return efficiency

def Effic_ULdigraph_Intermediate(N):
    """
    Smallest possible efficiency for strongly connected digraphs (ultra-long limit)
    containing exactly L = (N-1) + 1/2 N(N-1) arcs, or density = 1/2 + 1/N.

    Calculates the efficiency of a digraph consisting of the superposition of
    a directed ring with a complete directed acyclic graph whose arcs are
    pointing in the opposite orientations as the ring. That is, each vertex j
    points to all other vertices with i < j.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.

    Returns
    -------
    avpathlen : float
        The average pathlength of the digraph.

    See Also
    --------
    Effic_ULdigraph_Range1_MBS :
    Effic_ULdigraph_Range2 :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )

    # 1) CALCULATE THE PATHLENGTH
    Lmax = N*(N-1)
    term1 = N-1
    term2 = 0.5*N + EMconstant
    term3 = 0.0
    for i in range(1,N):
        term3 += scipy.special.psi(i+1)

    efficiency = ( term1*term2 + term3 ) / Lmax
    return efficiency

def Effic_ULdigraph_Disconnected(N,L):
    """
    Smallest possible efficiency of a digraph that is not strongly connected.

    Calculates the smallest possible efficiency that a digraph of arbitrary
    number of nodes N and number of directed arcs L could possibly have.
    When digraphs are allowed to be disconnected, the smallest efficiency any
    digraph could have is E_dUL = density. This is true for any value of
    L <= 1/2 N(N-1), or density <= 1/2. For larger link densities, special
    cases exists for which the smallest efficiency is larger than the density.
    While the configuration of such cases is as today not fully understood,
    the density is an excellent approximation for those cases as well.
    Therefore, this function returns, E_dUL = dens.

    Reference and citation
    ^^^^^^^^^^^^^^^^^^^^^^
    G. Zamora-López & R. Brasselet "Sizing complex networks" Comms Phys 2:144 (2019)

    Parameters
    ----------
    N : integer
        Number of nodes in the digraph.
    L : integer
        Number of undirected edges in the digraph.

    Returns
    -------
    efficiency : float
        The smallest efficiency a digraph could take.

    See Also
    --------
    Pathlen_ULdigraph :
    """
    # 0) SECURITY CHECKS
    if N < 2: raise ValueError( "Network needs at least two nodes, N > 1" )
    Ltot = N*(N-1)
    if L < 0:         raise ValueError( "L does not take negative values" )
    if L > int(Ltot): raise ValueError( "L out of range. max(L) = N*(N-1)" )

    # 1) CALCULATE THE EFFICIENCY
    efficiency = float(L) / Ltot
    return efficiency



##
