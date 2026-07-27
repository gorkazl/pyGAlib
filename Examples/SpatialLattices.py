# This file contains the functions developed for SiReNetA, to generate distance-
# aware network surrogates.
# THESE FUNCTIONS SHALL BECOME PART OF THE models.py MODULE AFTER TESTING,
# CLEANING AND VALIDATION.

"""

Spatially embedded (weighted) networks
--------------------------------------
SpatialWeightSorting
    Sorts the link weights of a network by the spatial distance between nodes.
SpatialLatticeFromNetwork
    Generates spatial weighted lattices with same weights as `con`.
"""

# Standard library imports

# Third party packages
import numpy as np
from numba import jit

# import galib
# from galib.models import*



## SPATIALLY EMBEDDED SURROGATES ###############################################
def SpatialWeightSorting(con, distmat, descending=True):
    """Sorts the link weights of a network by the spatial distance between nodes.

    The function reads the weights from a connectivity matrix and re-allocates
    them according to the euclidean distance between the nodes. The sorting
    conserves the position of the links, therefore, if `con` is a binary graph,
    the function will return a copy of `con`. The distance between nodes shall
    be given as input `distmat`.

    If descending = True, the larger weigths are assigned to the links between
    closer nodes, and the smaller weights to the links between distant nodes.

    If descending = False, the larger weights are assigned to the links between
    distant nodes, and the smaller weights to links between close nodes.

    Parameters
    ----------
    con : ndarray (2d) of shape (N,N).
        The connectivity matrix of the network.
    distmat : ndarray, rank-2.
        A matrix containing the spatial distance between all pair of ROIs.
        This can be either the euclidean distance, the fiber length or any
        other geometric distance.
    descending : boolean, optional.
        Determines whether links weights are assigend in descending or in
        ascending order, according to the euclidean distance between the nodes.

    Returns
    -------
    newcon : ndarray of rank-2 and shape (N x N).
        Connectivity matrix with weights sorted according to spatial distance
        between the nodes.

    """
    # 0) SECURITY CHECKS
    con_shape = np.shape(con)
    dist_shape = np.shape(distmat)
    if con_shape != dist_shape:
        raise ValueError( "Data not aligned. 'con' and 'distmat' of same shape expectted. " )

    # 1) EXTRACT THE NEEDED INFORMATION FROM THE con MATRIX
    N = len(con)
    # The indices of the links and their weights, distance
    nzidx = con.nonzero()
    weights = con[nzidx]
    distances = distmat[nzidx]

    # 2) SORT THE WEIGHTS IN DESCENDING ORDER
    weights.sort()
    if descending:
        weights = weights[::-1]

    # Get the indices that would sort the links by distance
    sortdistidx = distances.argsort()
    newidx = (nzidx[0][sortdistidx], nzidx[1][sortdistidx])

    # 3) CREATE THE NEW CONNECTIVITY WITH THE LINK WEIGHTS SORTED SPATIALLY
    newcon = np.zeros((N,N), np.float64)
    newcon[newidx] = weights

    return newcon

def SpatialLatticeFromNetwork(con, distmat, descending=True):
    """Generates spatial weighted lattices with same weights as `con`.

    The function reads the weights from a connectivity matrix and generates a
    spatially embedded weighted lattice, assigning the largest weights in
    descending order to the nodes that are closer from each other. Therefore,
    it requires also the euclidean distance between the nodes is given as input.

    If `con` is a binary graph of L links, the function returns a graph with
    links between the L spatially closest pairs of nodes.

    If `descending = True`, the larger weigths are assigned to the links between
    closer nodes, and the smaller weights to the links between distant nodes.

    If `descending = False`, the larger weights are assigned to the links between
    distant nodes, and the smaller weights to links between close nodes.

    Note
    ----
    Even if `con` is either a directed network or undirected but with asymmetric
    weights, the resulting lattice will be undirected and (quasi-)symmetric
    due to the fact that the spatial distance between two nodes is symmetric.

    Parameters
    ----------
    con : ndarray (2d) of shape (N,N).
        The connectivity matrix of the network.
    distmat : ndarray, rank-2.
        A matrix containing the spatial distance between all pair of ROIs.
        This can be either the euclidean distance, the fiber length or any
        other geometric distance.
    descending : boolean, optional.
        Determines whether links weights are assigend in descending or in
        ascending order, according to the euclidean distance between the nodes.

    Returns
    -------
    newcon : ndarray of rank-2 and shape (N x N).
        Connectivity matrix of a weighted lattice.

    """
    # 0) SECURITY CHECKS
    con_shape = np.shape(con)
    dist_shape = np.shape(distmat)
    if con_shape != dist_shape:
        raise ValueError( "Data not aligned. 'con' and 'distmat' of same shape expectted. " )

    # 1) EXTRACT THE NEEDED INFORMATION FROM THE con MATRIX
    N = len(con)

    # Sort the weights of the network
    weights = con.flatten()
    weights.sort()
    if descending:
        weights = weights[::-1]

    # Find the indices that sort the euclidean distances, from shorter to longer
    if descending:
        distmat[np.diag_indices(N)] = np.inf
    else:
        distmat[np.diag_indices(N)] = 0.0
    distances = distmat.ravel()
    sortdistidx = distances.argsort()
    newidx = np.unravel_index( sortdistidx, (N,N) )

    # And finally, create the coonectivity matrix with the weights sorted
    newcon = np.zeros((N,N), np.float64)
    newcon[newidx] = weights

    return newcon



##
