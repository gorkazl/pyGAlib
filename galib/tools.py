# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
MISCELLANOUS HELPER FUNCTIONS
=============================

This module contains miscellaneous helper functions useful for the analisys
of graphs and complex networks.

I/O AND DATA CONVERSIONS
------------------------
LoadFromPajek
    Reads a network from a text file with Pajek format.
Save2Pajek
    Saves a network into a text file with Pajek format.
LoadLabels
    Reads the labels of nodes from a text file.
SaveLabels
    Saves the labels of nodes into a text file
ReadPartition
    Reads a partition of nodes from a text file.
SavePartition
    Saves a partition of nodes into a text file.
ExtractSubmatrix
    Returns the sub-matrix composed by a set of nodes.
SymmetriseMatrix
    Converts a directed network into undirected by averaging the weights.
LaplacianMatrix
    Returns the Laplacian matrix of a given network.
CleanPaths
    Finds and removes in-place repeated, opposite paths from a list of paths.

ARRAY AND MATRIX COMPARISONS
----------------------------
ArrayCompare
    Compares whether two arrays are identical or not.
HammingDistance
    Computes the Hamming distance between two arrays of same shape.

ADDITIONAL MATH FUNCTIONS
-------------------------
NonZeroMin
    Returns the smallest non-zero value in an array.
CumulativeDistribution
    Computes the cumulative distribution of a dataset.
Factorial
    Computes the factorial of an integer number.
BinomialCoefficient
    Computes the binomial coefficient of n over m.
Quartiles
    Finds the Q1, Q2 and Q3 quartiles of a dataset.
AllPermutations
    Given a set, it returns all possible permutations.
AllCominations
    Given a set, finds all combinations of given size.
AllBipartitions
    Given a set, finds all its possible bipartitions.
MeanCorrelation
    Computes the Fisher-corrected mean value of correlation values.


...moduleauthor:: Gorka Zamora-Lopez <galib@zamora-lopez.xyz>

"""
from __future__ import division, print_function, absolute_import

import itertools
import types
import re
import warnings
import functools

import numpy as np

# from . import metrics
import galib.metrics


# I/O AND DATA CONVERSIONS ################################################
def LoadFromPajek(filepath, getlabels=False):
    """Reads a network from a text file with Pajek format.

    Parameters
    ----------
    filepath : string
        The source .net file of the network in Pajek format.
    getlabels : boolean, optional
        If True, the function also reads and returns the labels of the nodes.

    Returns
    -------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    labels : list, optional
        If getlabels = True, the function also return a list containing the
        names of the nodes.

    See Also
    --------
    LoadLabels, SaveLabels
    Save2Pajek : Saves a network into a file in Pajek-readable format.

    Notes
    -----
    1. The function automatically determines whether the network is directed,
    undirected and / or weighted.
    2. The returned adjacency matrix is of dtype 'int' or 'float', depending
    on the weights in the file.
    """
    # 0) OPEN THE FILE AND READ THE SIZE OF THE NETWORK
    pajekfile = open(filepath, 'r')
    firstline = pajekfile.readline()
    firstline = firstline.split()
    N = int(firstline[1])

    # 1) READ THE LABELS OF THE NODES IF WANTED
    if getlabels:
        labels = []

        # Security check, make sure that labels of nodes are listed in file
        line = pajekfile.readline()
        if line.split()[0] != '1':
            pajekfile.seek(1)
            print('LoadFromPajek() warning: No labels found to read.')

        # If labels are in file continue reading the labels.
        else:
            # If labels are wrapped in between quotes
            try:
                idx1 = line.index('"') + 1
                # Add the first label
                idx2 = line[idx1:].index('"')
                label = line[idx1:idx1+idx2]
                labels.append(label)

                # And now read the labels for the rest of the nodes
                for i in range(1,N):
                    line = pajekfile.readline()
                    idx1 = line.index('"') + 1
                    idx2 = line[idx1:].index('"')
                    label = line[idx1:idx1+idx2]
                    labels.append(label)

            # Otherwise, make a wild guess of what the label is
            except ValueError:
                # Add the first label
                label = line.split()[1]
                labels.append(label)

                # And now read the labels of the rest of the nodes
                for i in range(1,N):
                    line = pajekfile.readline()
                    label = line.split()[1]
                    labels.append(label)

    # 2) READ THE LINKS AND CREATE THE ADJACENCY MATRIX
    # 2.1) Find out whether the network is directed or undirected
    # while loop to skip empty lines if needed or the lines of the labels
    done = False
    while not done:
        line = pajekfile.readline()
        if line[0] == '*':
            if 'Edges' in line:
                directed = False
            elif 'Arcs' in line:
                directed = True
            else:
                print('Could not find whether network is directed or undirected')
                break
            done = True

    # 2.2) Read the first line contining a link
    line = pajekfile.readline()
    line = line.split()

    # If link information is BINARY, just read the adjacency list links
    if len(line) == 2:
        # 2.3) Declare the adjacency matrix and include the first link
        adjmatrix = np.zeros((N,N), np.uint8)
        i = int(line[0]) - 1
        j = int(line[1]) - 1
        adjmatrix[i,j] = 1
        if not directed:
            adjmatrix[j,i] = 1

        # 2.4) Include the rest of the links
        for line in pajekfile:
            i, j = line.split()
            i = int(i) - 1
            j = int(j) - 1
            adjmatrix[i, j] = 1
            if not directed:
                adjmatrix[j, i] = 1

    # If the link information is WEIGHTED, read the weighted links
    elif len(line) == 3:
        # 2.3) Find whether link weights are integer or floating poing
        i, j, aij = line
        outdtype = np.int
        try:
            outdtype(aij)
        except ValueError:
            outdtype = np.float

        # 2.4) Declare the adjacency matrix and include the first link
        adjmatrix = np.zeros((N, N), outdtype)
        i = int(i) - 1
        j = int(j) - 1
        adjmatrix[i, j] = outdtype(aij)
        if not directed:
            adjmatrix[j, i] = outdtype(aij)

        # 2.5) Read the rest of the file and fill-in the adjacency matrix
        for line in pajekfile:
            i, j, aij = line.split()
            i = int(i) - 1
            j = int(j) - 1
            adjmatrix[i, j] = outdtype(aij)
            if not directed:
                adjmatrix[j, i] = adjmatrix[i, j]

    # 3) CLOSE FILE AND RETURN RESULTS
    pajekfile.close()

    if getlabels:
        return adjmatrix, labels
    else:
        return adjmatrix

def Save2Pajek(filepath, adjmatrix, labels=[], directed=False, weighted=False):
    """Saves a network into a text file with Pajek format.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    filepath : string
        The target file in which the partition is to be saved.
    labels : list, optional
        A list containing the labels of each node.
    directed : Boolean, optional
        True if the network is directed, False otherwise.
    weighted : Boolean, optional
        If 'True', the weights of the links will be included in the file.
        Otherwise, it assumes the network is binary and it saves only the
        adjacency list, without information on the link weights.

    See Also
    --------
    LoadFromPajek : Loads a network from file in Pajek format.
    """
    # 0) SECURITY CHECK
    N = len(adjmatrix)
    if labels:
        if len(labels) != N:
            raise ValueError( "List of labels not aligned with network size" )

    # 1) OPEN THE TARGET FILE
    outfile = open(filepath, 'w')

    # 2) SAVE INFORMATION OF THE NODES
    print('*Vertices', N, file=outfile)

    # Wrie the list of nodes if labels have been given
    if labels:
        for i in range(N):
            # line = '%d "%s"' %(i+1, labels[i])
            line = '%d\t"%s"' % (i + 1, labels[i])
            print(line, file=outfile)

    # 3) SAVE THE LINKS
    # Save the links AND their WEIGHTS
    if weighted:
        # 3.1) Find whether weights are integers or floats
        if adjmatrix[0, 0].dtype in [np.uint8, np.uint, np.int8, np.int]:
            formatstring = '%d %d %d'
        elif adjmatrix[0, 0].dtype in [np.float16, np.float32, np.float, np.float64]:
            formatstring = '%d %d %f'

        # 3.2) Save the ARCS if directed
        if directed:
            print('*Arcs', file=outfile)
            for i in range(N):
                neighbours = adjmatrix[i].nonzero()[0]
                for j in neighbours:
                    line = formatstring % (i + 1, j + 1, adjmatrix[i, j])
                    print(line, file=outfile)

        # 3.2) Save the EDGES, if undirected
        else:
            print('*Edges', file=outfile)
            for i in range(N):
                neighbours = adjmatrix[i].nonzero()[0]
                for j in neighbours:
                    if j > i:
                        line = formatstring % (i + 1, j + 1, adjmatrix[i, j])
                        print(line, file=outfile)

    # Save ONLY the adjacency list
    else:
        formatstring = '%d %d'

        # 3.1) Save the ARCS if directed
        if directed:
            print('*Arcs', file=outfile)
            for i in range(N):
                neighbours = adjmatrix[i].nonzero()[0]
                for j in neighbours:
                    line = formatstring % (i + 1, j + 1)
                    print(line, file=outfile)

        # 3.1) Save the EDGES, if undirected
        else:
            print('*Edges', file=outfile)
            for i in range(N):
                neighbours = adjmatrix[i].nonzero()[0]
                for j in neighbours:
                    if j > i:
                        line = formatstring % (i + 1, j + 1)
                        print(line, file=outfile)

    # 4) CLOSE FILE AND FINISH
    outfile.close()

def LoadLabels(filepath):
    """Reads the labels of nodes from a text file.

    Parameters
    ----------
    filepath : string
        The path and filename in which the labels are stored.

    Returns
    -------
    data : list
        A list of strings with the label of each node.

    See Also
    --------
    SaveLabels : Saves the labels of nodes into a file
    """
    with open(filepath, 'r') as datafile:
        lines = [line.strip() for line in datafile.readlines()]
        # filter for empty lines
        return [line for line in lines if line]

def SaveLabels(filepath, labels):
    """Saves the labels of nodes into a text file.

    Parameters
    ----------
    labels : list-like
        A list or tuple of the labels of the nodes as strings.
    filepath : string
        The path and filename in which the labels will be stored.

    See Also
    --------
    LoadLabels : Reads the labels of nodes from a file
    """
    # 1) Create a string with all the text to be stored
    text = '\n'.join(labels)

    # 2) Open the datafile and save the text
    with open(filepath, 'w') as outfile:
        outfile.write(text)

def LoadPartition(filepath, sep=r'\s'):
    """Reads a partition of nodes from a text file.

    Every line of the file must contain the list of nodes in one community.

    Returns
    -------
    partition : list
        A list containing as many lists of nodes as communities.

    See Also
    --------
    SavePartition : Saves a partition of nodes into a file.
    """

    # 0) READ THE DATA FILE
    with open(filepath, 'r') as datafile:
        partition = []
        for line in datafile:
            # 1) CUT AT SEPARATOR
            newcom = re.split(sep, line.strip())
            # 2) CONVERT TO DECIMALS AND APPEND TO LIST
            partition.append([int(x) for x in newcom])

        return partition

def SavePartition(filepath, partition, sep=' '):
    """Saves a partition of nodes into a text file.

    Parameters
    ----------
    partition : list
        A list containing as many lists of nodes as communities.
    filepath : string
        The target file in which the partition is to be saved.
    labels : boolean, optional

    See Also
    --------
    LoadPartition : Reads a partition of nodes saved in a text file.
    """
    def _com_to_string(com):
        """ helper function to string format node weights """
        return sep.join([str(x) for x in com])

    # joining lines with line-breaks
    text = '\n'.join([_com_to_string(com) for com in partition])

    # writing formatted matrix to file
    with open(filepath, 'w') as outfile:
        outfile.write(text)

def ExtractSubmatrix(adjmatrix, nodelist1, nodelist2=[]):
    """Returns the sub-matrix composed by a set of nodes.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2.
        The matrix from which the submatrix will be extracted.
    nodelist1 : arra_like. list, tuple, set or ndarray of rank-1.
        A list of indices with the rows of the matrix to be extracted.
        Prefered type is 'ndarray'. All other sequence types (lists, tuples and
        sets) will be converted to an ndarray of dtype = numpy.int.
    nodelist2 : arra_like. list, tuple, set or ndarray of rank-1 (optional)
        If not given (nodelist2=None), then nodelist2 is nodelist1.
        A list of indices with the columns of the matrix to be extracted.
        Prefered type is 'ndarray'. All other sequence types (lists, tuples and
        sets) will be converted to an ndarray of dtype = numpy.int.

    Returns
    -------
    A submatrix of 'adjmatrix' of size n1 x n2, where n1 and n2 are the
    lengths of nodelist1 and nodelist2 respectively.

    Examples
    --------
    >>> testmat = arange(1,17).reshape(4,4)
    >>> ExtractSubmatrix(testmat,[0,3])
    array([[ 1,  4],
           [13, 16]])
    >>>
    >>> ExtractSubmatrix(testmat,[0,3],[1,2,3])
    array([[ 2,  3,  4],
           [14, 15, 16]])

    """
    # 0) CHECK WHETHER LISTS OF NODES ARE GIVEN AS ARRAYS. OTHERWISE, CONVERT.
    # Check nodelist1
    if type(nodelist1) != np.ndarray:
        # warnings.simplefilter('once', UserWarning)
        # warnings.warn("Prefered type for parameter 'nodelist1' is numpy.ndarray. Lists, tuples and sets are allowed but are converted to ndarrays.", \
        #             stacklevel=2)

        if type(nodelist1) == set:
            nodelist1 = list(nodelist1)
        if type(nodelist1) in [list, tuple]:
            nodelist1 = np.array(nodelist1,np.int)

    # Check nodelist2
    if len(nodelist2) == 0:
        nodelist2 = nodelist1.copy()
    else:
        if type(nodelist2) != np.ndarray:
            # warnings.simplefilter('once', UserWarning)
            # warnings.warn("Prefered type for parameter 'nodelist2' is numpy.ndarray. Lists, tuples and sets are allowed but are converted to ndarrays.", \
            #                 stacklevel=2)

            if type(nodelist2) == set:
                nodelist2 = list(nodelist2)
            if type(nodelist2) in [list, tuple]:
                nodelist2 = np.array(nodelist2,np.int)

    # 1) CREATE LISTS OF INDICES FOR SLICING
    N1 = len(nodelist1)
    N2 = len(nodelist2)
    xindices = np.zeros(N1*N2, np.int)
    for ncounter, node in enumerate(nodelist1):
        startidx = ncounter*N2
        endidx = startidx + N2
        xindices[startidx:endidx] = node

    yindices = np.zeros(N1*N2, np.int)
    for n in range(N1):
        startidx = n*N2
        endidx = startidx + N2
        yindices[startidx:endidx] = nodelist2

    # 2) EXTRACT THE SUBMATRIX AND FINISH
    return adjmatrix[xindices, yindices].reshape(N1, N2)

def SymmetriseMatrix(adjmatrix):
    """Converts a directed network into undirected by averaging the weights,
    i.e., bij = 1/2 * (aij + aji).

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.

    Returns
    -------
    An adjacency matrix of the same shape, of dype=float, with values
    """

    if galib.metrics.Reciprocity(adjmatrix) == 1:
    # if Reciprocity(adjmatrix) == 1:
        return adjmatrix
    else:
        return 0.5 * (adjmatrix + adjmatrix.T)

def LaplacianMatrix(adjmatrix):
    """Returns the Laplacian matrix of a given network.

    The Laplacian matrix of a network is the matrix L = D - A, where A is
    the adjacency matrix and D is a diagonal matrix with elements Dii being
    the corresponding row sum of A.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.

    Returns
    -------
    laplacianmatrix : ndarray of same shape and dtype as adjmatrix

    Warning
    -------
    If the network is directed, L should be computed out of the transposed
    adjacency matrix, otherwise it would represent the dynamics of the nodes
    driven by their outputs, not by their inputs. We recommend to introduce
    adjmatrix.T into the function.
    """
    if adjmatrix.dtype in [np.uint, np.uint0, np.uint8, np.uint16, np.uint32, np.uint64]:
        adjmatrix = adjmatrix.astype(int)
    N = len(adjmatrix)

    laplacianmatrix = np.identity(N, dtype=adjmatrix.dtype) * adjmatrix.sum(axis=1)
    laplacianmatrix -= adjmatrix

    return laplacianmatrix

def CleanPaths(pathlist):
    """Finds and removes in-place repeated, opposite paths from a list of paths.

    Parameters
    ----------
    pathlist : list
        A list containing paths as list or array-like.

    See Also
    --------
    PathsAllinOne : Computes several graph distance and path measures at once.
    CleanCycles : Removes repeated cycles from a list of cycles.

    Notes
    -----
    1. The function modifies the parameter 'pathlist' in-place and it does
    not return a value. On the contrary, the function CleanCycles() does
    return another list.
    2. The function is intended to be used with the output of PathsAllinOne(),
    and ONLY in the case of undirected graphs because PathsAllinOne() returns
    all paths twice, from i to j and from j to i. Do not apply CleanPaths()
    to the output from directed networks.

    Examples
    --------
    >>> paths = [[1,2,3],[3,2,1],[2,4,7,8],[8,7,4,2]]
    >>> CleanPaths(paths)
    >>> paths
    [[1, 2, 3], [2, 4, 7, 8]]
    """
    for path1 in pathlist:
        for path2 in pathlist[::-1]:
            if path2[::-1] == path1:
                pathlist.remove(path2)
                break

## ARRAY AND MATRIX COMPARISONS ###############################################
def HammingDistance(array1, array2, normed=False):
    """Computes the Hamming distance between two arrays of same shape.

    The Hamming distance counts the fraction of elements that are different
    in two arrays. It is 0 when all the two arrays are identical and 1 when
    the two arrays have no common elements.

    Note
    ----
    In case of binary data, it can be normed to account for the random
    probability of coincidences in two arrays with n_1,0 and n_2,0 number
    of zero entries and n_1,1 and n_2,1 number of ones.

    Parameters
    ----------
    array1 : ndarray of any shape and dimension
    array2 : ndarray of same shape and dimension as array1
    normed : Boolean, optional
        If True, the function computes the expected Hamming distance between
        the two vectors. The expected value is based on the expected number
        of coincidences that would be found if both arrays were random
        binary vectors with n_1 and n_2 number of ones respectively (giving
        length - n_1 and length - n_2 number of zeros).

    Returns
    -------
    h_dist : floating point.
        The Hamming distance between the two arrays.
    exp_hdist : floating point. Optional, only if 'expected = True'
        The expected Hamming distance between the vectors
    """
    # 0) PREPARE FOR CALCULATIONS
    # 0.1) Convert the arrays into rank-1 arrays
    if len(np.shape(array1)) > 1:
        array1 = array1.reshape(-1)
    if len(np.shape(array2)) > 1:
        array2 = array2.reshape(-1)

    # 0.2) Security check
    if len(array1) != len(array2):
        raise ValueError( "Arrays are not aligned" )

    # 1) COUNT THE NUMBER OF COINCIDENCES
    similarity = (array1 == array2)
    n_equal = similarity.sum()

    # 2) COMPUTE THE HAMMING DISTANCE
    length = len(array1)
    h_dist = 1. - float(n_equal) / length

    # 3) RETURN RESULT ACCORDING TO OPTIONS
    # Standard Hamming distance
    if not normed:
        return h_dist

    # Normalized Hamming distance
    else:
        # Count the number of ones in the two arrays
        n_1 = len(array1.nonzero()[0])
        n_2 = len(array2.nonzero()[0])

        # Estimate the expected number of random coincidences
        exp_nc = 1.0 / length * (n_1 * n_2 + (length - n_1) * (length - n_2))

        # The expected Hamming distance
        exp_hdist = 1.0 - exp_nc / float(length)

        return h_dist, exp_hdist

## SOME ADDITIONAL MATH ######################################################
def NonZeroMin(data):
    """Returns the smallest non-zero value in an array.

    Parameters
    ----------
    data : array-like
        A list, tuple or array of numbers.

    Returns
    -------
    An integer or real value, depending on data's dtype.
    """

    # 1) Convert lists and tuples into arrays
    if type(data) != 'numpy.ndarray':
        data = np.array(data)

    # 2) Find the minimal non-zero value and return.
    idx = np.where(data)
    return data[idx].min()

def CumulativeDistribution(data, nbins, range=None, normed=True, centerbins=False):
    """Computes the cumulative distribution of a dataset.

    Parameters
    ----------
    data : array_like
        Input data. The histogram is computed over the flattened array.
    nbins : int or float
        If `nbins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
    range : (float, float), optional
        The lower and upper range of the bins. If not provided, range
        is simply (data.min(), data.max()). Values outside the range are
        ignored.
    normed : bool, optional
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that the sum of the
        histogram values will not be equal to 1 unless bins of unity
        width are chosen; it is not a probability *mass* function.
    centerbins : bool, optional
        If False, the positions of the bins ('xdata' variable) returned
        correspond to the lower limit of the each bin, as it is returned
        originally by the histogram() function of numpy. If True, 'xdata' will
        contain the value at the center of the bin.

    Returns
    -------
    xdata : ndarray
        Returns the bin location. This will be the left edge of the bins if
        option 'centerbins=False', and the center of the bin if option
        centerbins=True.
    ydata : ndarray
        The values of the cumulative probability.

    See Also
    --------
    ''histogram()'' function for further documentation.
    """

    # 1) COMPUTE THE DISTRIBUTION OF THE DATA
    ydata, xdata = np.histogram(data, nbins, range, normed)

    # 1.1) Compute the cumulative sum of the probability
    ydata = ydata.cumsum()

    # 2) RETURN THE RESULTS
    if centerbins:
        dif = 0.5 * (xdata[-1] - xdata[0]) / nbins
        xdata += dif

    if normed:
        norm = 1.0 / ydata[-1]
        ydata *= norm

        return xdata[:-1], ydata

    else:
        return xdata[:-1], ydata

def Factorial(x):
    """Computes the factorial of an integer number.
    """
    # 0) SECURITY CHECK
    if not isinstance(x, int):
        raise ValueError( "'Factorial' function only accepts integers" )

    # 1) COMPUTE THE FACTORIAL
    if x == 0 or x == 1:
        return 1
    else:
        return functools.reduce(lambda x, y: x * y, range(1, x + 1))

def BinomialCoefficient(n, m):
    """Computes the binomial coefficient of n over m.
    """
    if m == 0:
        return 1

    elif m == 1:
        return n

    else:
        ma = max(n - m, m)
        mi = min(n - m, m)

        enum = functools.reduce(lambda x, y: x * y, range(ma + 1, n + 1), 1)

        return enum / Factorial(mi)

def Quartiles(data):
    """
    Finds the 1st, 2nd and 3rd quartiles of a dataset.

    Parameters
    ----------
    data : array_like
        Input data.

    Returns
    -------
    Q1 : float. The first quartile.
    Q2 : float. The second quartile or median.
    Q3 : float. The third quartile.

    See Also
    ---------
    numpy.percentile
    """
    q = np.percentile(data, [25, 50, 75])

    return q[0], q[1], q[2]

def AllPermutations(data):
    """Given a set, it returns all possible permutations.

    Parameters
    ----------
    data : array_like
        A list, tuple or array containing data.

    Returns
    -------
    permutationlist : list
        A list containing all possible permutations of the elements in data.

    See Also
    --------
    itertools.permutations, AllCombinations, AllBipartitions

    Examples
    --------
    >>> AllPermutations((1,2,3))
    >>> [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]

    >>> AllPermutation(('test','case',3))
    >>> [('test', 'case', 3),
        ('test', 3, 'case'),
        ('case', 'test', 3),
        ('case', 3, 'test'),
        (3, 'test', 'case'),
        (3, 'case', 'test')]
    """
    if len(data) <= 1:
        return data

    return [p for p in itertools.permutations(data)]

def AllCombinations(data, comblength):
    """Given a set, finds all combinations of given size.

    Returns a list with all possible permutations of length <= comblength
    of the elements in data. If N is the length of the data, the number
    of combinations of length 'comblength' is the binomial coefficient
    BinomialCoefficient(N,comblength).

    Parameters
    ----------
    data : array_like
        A list, tuple or array containing data.
    comblength : int
        Number of elements of the permutation lists.

    Returns
    -------
    combinationlist : list
        A list of all possible permutations of length 'comblength' with the
        elements in 'data'.

    See Also
    --------
    itertools.combinations, AllPermutations, AllBipartitions

    Examples
    --------
    >>>
    >>> [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    >>>
    >>> AllCombinations((1,2,3,4), 3)
    >>> [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    >>>
    >>> AllCombinations(('s','p','a','m'), 3)
    >>> [('s', 'p', 'a'), ('s', 'p', 'm'), ('s', 'a', 'm'), ('p', 'a', 'm')]
    """
    return [c for c in itertools.combinations(data, comblength)]

def AllBipartitions(data, comblength=None):
    """Given a set, finds all possilbe bipartitions.

    If optional 'comblength=None', it will return all bipartitions of
    all sizes. If 'comblength' is an integer, it will return all bipartitions
    in which one of the subsets has size N1 = comblength and complementary
    has size N2 = len(data) - comblength.

    Parameters
    ----------
    data : array_like
        A list, tuple or array containing data.
    comblength : integer, optional
        The size of one of the two complementary subsets.

    Returns
    -------
    result : list
        A list with all possible bipartitions of the elements in 'data'.

    See Also
    --------
    AllCombinations, AllPermutations

    Examples
    --------
    >>> AllBipartitions((1,2,3,4))
    >>> [((1,), (2, 3, 4)),
        ((2,), (1, 3, 4)),
        ((3,), (1, 2, 4)),
        ((4,), (1, 2, 3)),
        ((1, 2), (3, 4)),
        ((1, 3), (2, 4)),
        ((1, 4), (2, 3))]
    >>>
    >>> AllBipartitions(('s','p','a','m'))
    >>> [(('s',), ('a', 'p', 'm')),
        (('p',), ('a', 's', 'm')),
        (('a',), ('p', 's', 'm')),
        (('m',), ('a', 'p', 's')),
        (('s', 'p'), ('a', 'm')),
        (('s', 'a'), ('p', 'm')),
        (('s', 'm'), ('a', 'p'))]
    """
    # 0) PREPARE FOR CALCULATIONS
    Ndata = len(data)
    iseven = Ndata % 2 == 0
    Nmax = Ndata // 2
    setdata = set(data)

    # Security check, avoid confusing comblength=0 with None
    if comblength == 0:
        raise ValueError( "'comblength' out of range, please make sure 0 < comblength < len(data)." )

    ########################################################################
    # 1) DEFINE THE TWO FUNCTIONS THAT WILL DO THE ACTUAL JOB
    def AllBipartitionsAllSizes(data):
        """Given a dataset, it finds all bipartitions of all sizes.
        """
        # 1) FIND ALL THE COMBINATIONS UP TO LENGTH < len(data)/2
        bipartitions = []
        for n in range(1, Nmax):
            # 1.1) Find all combinations of given size out of dataset
            combinations = [c for c in itertools.combinations(data, n)]
            # 1.2) Sort and find complementary sets
            for comb in combinations:
                complementary = tuple(setdata - set(comb))
                bipartitions.append((comb, complementary))

        # 2) FIND AND SORT THE BIPARTITIONS OF SIZE Nmax
        # 2.1) Find all combinations of size Nmax
        combinations = [c for c in itertools.combinations(data, Nmax)]
        # 2.2) Sort and find complementary sets
        ncombs = len(combinations)
        # Ignore repeated combinations if both subsets are of size = Nmax
        if iseven: ncombs = ncombs // 2

        for i in range(ncombs):
            comb = combinations[i]
            combset = set(comb)
            complementary = setdata - combset
            bipartitions.append((comb, tuple(complementary)))

        return bipartitions

    def AllBipartitionsGivenSize(data, comblength):
        """Given a dataset, finds all bipartitions in which one of the
        subsets is of size 'comblength'. The size of the complementary
        subset is N - 'comblength', where N is the size of the data.
        """
        # 0) SECURITY CHECKS
        comblength = int(comblength)
        if (comblength < 1 or comblength >= Ndata):
            raise ValueError( "'comblength' out of range, please make sure 0 < comblength < len(data)." )

        # 1) FIND ALL THE COMBINATIONS OF SIZE 'comblength'
        combinations = [c for c in itertools.combinations(data, comblength)]

        # 2) SORT AND FIND THE COMPLEMENTARY SETS
        bipartitions = []

        # If size of set is even, and we are looking for bipartitions of size
        # Nmax, then we need to remove repeated combinations
        if iseven and comblength == Nmax:
            ncombs = len(combinations)
            for i in range(ncombs // 2):
                comb = combinations[i]
                combset = set(comb)
                complementary = setdata - combset
                bipartitions.append((comb, tuple(complementary)))

        # For any other case, bipartitions are found as usual
        else:
            for comb in combinations:
                complementary = tuple(setdata - set(comb))
                bipartitions.append((comb, complementary))

        return bipartitions

    ########################################################################
    # 2) DO THE JOB
    if comblength:
        result = AllBipartitionsGivenSize(data, comblength)
    else:
        result = AllBipartitionsAllSizes(data)

    return result

def MeanCorrelation(data, tolerance=10**(-15)):
    """Computes the Fisher-corrected mean value of correlation values.

    Parameters
    ----------
    data : array-like ndarray or list
        A tuple, list or array containing the values of correlation that
        are to be averaged. Accepts only 1-dimensional data, no matrices.
    tolerance : float
        Small jitter in the correlation values around exact 1.0 and -1.0
        that will also be accepted. Largest correlation values are
        set to 1.0 - tolerance and smallest correlation values are set to
        -1.0 + tolerance. This avoids +/- inf values to be returned by the
        arctanh() function. A single inf value would totally bias the mean.

    Returns
    -------
    float value
        Mean value of correlation in the range (-1,1).
    """
    # Security checks and data preparation
    dataset = np.array(data, np.float64)
    if len(np.shape(dataset)) > 1:
        dataset = dataset.flatten()

    if dataset.max() > 1.0 + tolerance:
        raise ValueError( "Correlation values larger than 1.0 in data" )
    if dataset.min() < -1.0 - tolerance:
        raise ValueError( "Correlation values smaller than -1.0 in data" )

    # Remove extremal values to avoid infinite values in the arctanh(x).
    idx = np.where(data >= 1.0)
    if idx:
        data[idx] = 1.0 - tolerance

    idx = np.where(data <= -1.0)
    if idx:
        data[idx] = -1.0 + tolerance

    # Compute the Fisher corrected
    newdata = np.arctanh(data)
    newdata = np.where(newdata==np.inf, 10,newdata)
    newdata = np.where(newdata==-np.inf, -10,newdata)

    return np.tanh(newdata.mean())

#
