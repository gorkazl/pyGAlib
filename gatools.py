"""
=============================
MISCELLANOUS HELPER FUNCTIONS
=============================

This module contains miscellaneous helper functions useful for the analisys
of graphs and complex networks.

I/O AND DATA CONVERSIONS
========================
LoadLabels
    Reads the labels of nodes from a text file.
SaveLabels
    Saves the labels of nodes into a text file
ReadPartition
    Reads a partition of nodes from a text file.
SavePartition
    Saves a partition of nodes into a text file.
LoadFromPajek
    Reads a network from a text file with Pajek format.
Save2Pajek
    Saves a network into a text file with Pajek format.
ExtractSubmatrix
    Returns the sub-matrix composed by a set of nodes.
SymmetriseMatrix
    Converts a directed network into undirected by averaging the weights.
LaplacianMatrix
    Returns the Laplacian matrix of a given network.
CleanPaths
    Finds and removes in-place repeated, opposite paths from a list of paths.

ARRAY AND MATRIX COMPARISONS
============================
ArrayCompare
    Compares whether two arrays are identical or not.
HammingDistance
    Computes the Hamming distance between two arrays of same shape.

ADDITIONAL MATH FUNCTIONS
=========================
NonZeroMin
    Returns the smallest non-zero value in an array.
CumulativeDistribution
    Computes the cumulative distribution of a dataset.
Factorial
    Computes the factorial of an integer number.
BinomialCoefficient
    Computes the binomial coefficient of n over m.
StdDeviation
    Returns the mean value and standard deviation of a dataset.
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
"""

__author__ = "Gorka Zamora-Lopez"
__email__ = "Gorka.zamora@ymail.com"
__copyright__ = "Copyright 2013-2015"
__license__ = "GPL"
__update__ = "07/01/2015"

import itertools
import re

import galib
import numpy as np


## I/O AND DATA CONVERSIONS ################################################
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
            partition.append([float(x) for x in newcom])

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
    1. The function automatically determines whether the network is directed
    or undirected.
    2. The returned adjacency matrix is of dtype 'int' or 'float', depending
    on the weights in the file. The actual meaning of 'int' and 'float'
    depends on the machine. In 64-bit machines they will be 'int' = 'int64'
    and 'float' = 'float64'. See the NumPy documentation for further details.
    For memory savings in the case of binary networks the returned adjacency
    matrix shall be converted to 'uint8' using array.astype(uint8) method.
    The function automatically determines whether the network is directed
    or undirected.
    """
    # 0) OPEN THE FILE AND READ THE SIZE OF THE NETWORK
    pajekfile = open(filepath, 'r')
    firstline = pajekfile.readline()
    firstline = firstline.split()
    N = int(firstline[1])

    # 1) READ THE LABELS OF THE NODES IF SO DESIRED
    if getlabels:
        labels = []
        for i in xrange(N):
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
                print 'Could not assert whether network is directed or undirected'
                break
            done = True

    # 2.2) Find whether links are real valued or integer
    line = pajekfile.readline()
    i, j, aij = line.split()

    outdtype = np.int
    try:
        outdtype(aij)
    except ValueError:
        outdtype = np.float

    # 2.3) Declare the adjacency matrix and include the first link
    adjmatrix = np.zeros((N, N), outdtype)
    i = int(i) - 1
    j = int(j) - 1
    adjmatrix[i, j] = outdtype(aij)
    if not directed:
        adjmatrix[j, i] = outdtype(aij)

    # 2.4) Read the rest of the file and fill-in the adjacency matrix
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

def Save2Pajek(filepath, adjmatrix, labels=[], directed=True):
    """Saves a network into a text file with Pajek format.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network.
    filepath : string
        The target file in which the partition is to be saved.
    labels : list, optional
        A list containing the labels of each node.
    directed: Boolean, optional
        True if the network is directed, False otherwise.

    See Also
    --------
    LoadFromPajek : Loads a network from file in Pajek format.
    """
    # 0) SECURITY CHECK
    N = len(adjmatrix)
    if labels:
        assert len(labels) == N, 'List of labels not aligned with network size'

    # 1) OPEN THE TARGET FILE
    outfile = open(filepath, 'w')

    # 2) SAVE INFORMATION OF THE NODES
    print >> outfile, '*Vertices', N

    # 2.1) If the labels of the nodes have been given
    if labels:
        for i in xrange(N):
            # line = '%d "%s"' %(i+1, labels[i])
            line = '%d %s' % (i + 1, labels[i])
            print >> outfile, line

    # 2.2) If the labels of the nodes are not given
    else:
        for i in xrange(N):
            line = '%d "%d"' % (i + 1, i + 1)
            print >> outfile, line

    # 3) SAVE THE LINKS
    if adjmatrix[0, 0].dtype in [np.uint8, np.uint, np.int8, np.int]:
        formatstring = '%d %d %d'
    elif adjmatrix[0, 0].dtype in [np.float32, np.float, np.float64]:
        formatstring = '%d %d %f'

    # 3.1) Save the ARCS if directed
    if directed:
        print >> outfile, '*Arcs'
        for i in xrange(N):
            neighbours = adjmatrix[i].nonzero()[0]
            for j in neighbours:
                line = formatstring % (i + 1, j + 1, adjmatrix[i, j])
                print >> outfile, line

        # Close the outfile and finish
        outfile.close()

    # 3.2) Save the EDGES, if undirected
    else:
        print >> outfile, '*Edges'
        for i in xrange(N):
            neighbours = adjmatrix[i].nonzero()[0]
            for j in neighbours:
                if j > i:
                    line = formatstring % (i + 1, j + 1, adjmatrix[i, j])
                    print >> outfile, line

        # Close the file and finish
        outfile.close()

def ExtractSubmatrix(adjmatrix, nodelist1, nodelist2=None):
    """Returns the sub-matrix composed by a set of nodes.

    Parameters
    ----------
    adjmatrix : ndarray of rank-2.
    nodelist1 : list, tuple or ndarray.
        A list of indices meaning the rows to be extracted.
    nodelist2 : list, tuple or ndarray. Optional.
        A list of indices meaning the columns to be extracted.
        If a list not given (nodelist2=None), then nodelist2 is nodelist1

    Returns
    -------
    A submatrix of 'adjmatrix' of size n1 x n2, where n1 and n2 are the
    lengths of nodelist1 and nodelist2 respectively.

    Examples
    --------
    >>> testmat = array(((1,2,3),(4,5,6),(7,8,9)),int)
    >>> ExtractSubmatrix(testmat,[0,3])
    array([[ 1,  4],
           [13, 16]])
    >>>
    >>> ExtractSubmatrix(testmat,[0,3],[1,2,3])
    array([[ 2,  3,  4],
           [14, 15, 16]])

    """
    # 0) CHECK WHETHER LISTS OF NODES ARE GIVEN AS ARRAYS
    if type(nodelist1) == np.ndarray:
        nodelist1 = list(nodelist1)
    if nodelist2 == None:
        nodelist2 = nodelist1
    else:
        if type(nodelist2) == np.ndarray:
            nodelist2 = list(nodelist2)

    # 1) CREATE LISTS OF INDICES FOR SLICING
    N1 = len(nodelist1)
    N2 = len(nodelist2)
    xindices = []
    for node in nodelist1:
        xindices += [node] * N2
    yindices = nodelist2 * N1

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

    if galib.Reciprocity(adjmatrix) == 1:
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
    length = len(adjmatrix)

    laplacianmatrix = np.identity(length, dtype=adjmatrix.dtype) * adjmatrix.sum(axis=1)
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
    flat_array1 = array1.flatten()
    flat_array2 = array2.flatten()

    # 0.2) Security check
    assert len(flat_array1) == len(flat_array2), 'Arrays are not aligned'

    # 1) COUNT THE NUMBER OF COINCIDENCES
    similarity = (flat_array1 == flat_array2)
    n_equal = similarity.sum()

    # 2) COMPUTE THE HAMMING DISTANCE
    h_dist = 1. - float(n_equal) / len(flat_array1)

    # 3) RETURN RESULT ACCORDING TO OPTIONS
    # Standard Hamming distance
    if not normed:
        return h_dist

    # Normalized Hamming distance
    else:
        # Count the number of ones in the two arrays
        n_1 = len(flat_array1.nonzero()[0])
        n_2 = len(flat_array2.nonzero()[0])

        # Estimate the expected number of random coincidences
        length = len(flat_array1)
        exp_nc = 1. / length * (n_1 * n_2 + (length - n_1) * (length - n_2))

        # The expected Hamming distance
        exp_hdist = 1. - exp_nc / float(length)

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


def CumulativeDistribution(data, nbins, range=None, normed=True):
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

    Returns
    -------
    ydata : ndarray
        The values of the cumulative probability.
    xdata : ndarray
        Returns the bin edges.

    See Also
    --------
    ''histogram()'' function for further documentation.
    """

    # 1) COMPUTE THE DISTRIBUTION OF THE DATA
    ydata, xdata = np.histogram(data, nbins, range, normed)

    # 1.1) Compute the cumulative sum of the probability
    ydata = ydata.cumsum()

    # 2) RETURN THE RESULTS
    if normed:
        norm = 1. / ydata[-1]
        ydata *= norm

        return xdata[:-1], ydata

    else:
        return xdata[:-1], ydata

def Factorial(x):
    """Computes the factorial of an integer number.
    """
    # 0) SECURITY CHECK
    assert isinstance(x, int), "'factorial' function only accepts integers"

    # 1) COMPUTE THE FACTORIAL
    if x == 0 or x == 1:
        return 1

    else:
        return reduce(lambda x, y: x * y, xrange(1, x + 1))

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

        enum = reduce(lambda x, y: x * y, xrange(ma + 1, n + 1), 1)

        return enum / Factorial(mi)

def StdDeviation(data):
    """Returns the mean value and the standard deviation of an array of data.
    It is a simple wrapper using the numpy ufuncs a.mean() and a.dev().
    """
    if type(data) == np.ndarray:
        return data.mean(), data.std()
    else:
        data = np.array(data)
        return data.mean(), data.std()

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
    >>> AllCombinations((1,2,3,4), 2)
    >>> [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    >>>
    >>> AllCombinations((1,2,3,4), 3)
    >>> [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    >>>
    >>> AllCombinations(('s','p','a','m'), 3)
    >>> [('s', 'p', 'a'), ('s', 'p', 'm'), ('s', 'a', 'm'), ('p', 'a', 'm')]
    """
    return [c for c in itertools.combinations(data, comblength)]

def AllBipartitions(data):
    """Given a set, finds all its possilbe bipartitions.

    Parameters
    ----------
    data : array_like
        A list, tuple or array containing data.

    Returns
    -------
    bipartitions : list
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
    # 0) FIND WHETHER LENGTH OF DATASET IS ODD OR EVEN
    Ndata = len(data)
    iseven = Ndata % 2 == 0

    # 1) FIND ALL THE COMBINATIONS UP TO LENGTH len(data)/2
    Nmax = Ndata / 2
    setdata = set(data)
    bipartitions = []
    for n in xrange(1, Nmax):
        combinations = AllCombinations(data, n)
        for comb in combinations:
            complementary = tuple(setdata - set(comb))
            bipartitions.append((comb, complementary))

    # 2) ORGANIZSE THE COMBINATIONS INTO BIPARTITIONS
    combinations = AllCombinations(data, Nmax)

    helperlist = []
    for comb in combinations:
        combinationset = set(comb)
        if iseven and combinationset in helperlist:
            continue
        else:
            complementary = tuple(setdata - combinationset)
            bipartitions.append((comb, complementary))
            helperlist.append(set(complementary))

    return bipartitions

def MeanCorrelation(data, tolerance=10 ** (-15)):
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
    assert len(np.shape(data)) == 1, 'Only 1D arrays or lists of data accepted'
    dataset = np.array(data, np.float64)

    assert dataset.max() < 1.0 + tolerance, 'Correlation values larger than 1.0 in data'
    assert dataset.min() > -1.0 - tolerance, 'Correlation values smaller than -1.0 in data'

    # Remove extremal values to avoid infinite values in the arctanh(x).
    idx = np.where(data >= 1.0)
    if idx:
        data[idx] = 1.0 - tolerance

    idx = np.where(data <= -1.0)
    if idx:
        data[idx] = -1.0 + tolerance

    # Compute the Fisher corrected
    newdata = np.arctanh(data)

    return np.tanh(newdata.mean())
