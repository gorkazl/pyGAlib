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
"""

import itertools
from numpy import*
from galib import Reciprocity


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
    # 1) Read the data file
    datafile = open(filepath, 'r')
    data = datafile.read()
    datafile.close()
    
    # 2) Arange the data into a list
    data = data.split('\n')
    
    # 3) Check if last line is empty
    if not data[-1]:
        del data[-1]
    
    return data
    
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
    text = ''
    for node in labels:
        text += (node + '\n')

    # 2) Open the datafile and save the text
    outfile = open(filepath, 'w')
    outfile.write(text)
    outfile.close()
    
    # 3) Clean some trash
    del text

def LoadPartition(filepath):
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
    datafile = open(filepath, 'r')
    lines = datafile.readlines()
    datafile.close()
    Ncoms = len(lines)
    
    # 1) ORGANIZE THE DATA INTO LISTS
    partition = []
    for i in xrange(Ncoms):
        newcom = lines[i].split()
        newcom = array(newcom, int).tolist()
        partition.append(newcom)

    return partition

def SavePartition(filepath, partition):
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
    # 0) Open the target file
    outfile = open(filepath, 'w')
    
    # 1) Organize the data of each community into a text string
    for com in partition:
        line = ''
        for node in com:
            line += '%d ' %node
        # 2) Save the string to file
        print >> outfile, line

    outfile.close()  

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
    1. The function authomatically determines whether the network is directed
    or undirected.
    2. The returned adjacency matrix is of dtype 'int' or 'float', depending
    on the weights in the file. The actual meaning of 'int' and 'float'
    depends on the machine. In 64-bit machines they will be 'int' = 'int64'
    and 'float' = 'float64'. See the NumPy documentation for further details.
    For memory savings in the case of binary networks the returned adjacency
    matrix shall be converted to 'uint8' using array.astype(uint8) method.
    The function authomatically determines whether the network is directed
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
            #lowindx = line.index('"') + 1
            #highindx = line[lowindx:].index('"')
            #name = line[lowindx:lowindx+highindx]
            #labels.append(name)
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
    
    outdtype = int
    try:
        outdtype(aij)
    except ValueError:
        outdtype = float

    # 2.3) Declare the adjacency matrix and include the first link
    adjmatrix = zeros((N,N),outdtype)
    i = int(i) - 1
    j = int(j) - 1
    adjmatrix[i,j] = outdtype(aij)
    adjmatrix[j,i] = outdtype(aij)    
    
    # 2.4) Read the rest of the file and fill-in the adjacency matrix
    for line in pajekfile:
        i, j, aij = line.split()
        i = int(i) - 1
        j = int(j) - 1
        adjmatrix[i,j] =  outdtype(aij)
        if not directed:
            adjmatrix[j,i] = adjmatrix[i,j]
    
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
            #line = '%d "%s"' %(i+1, labels[i])
            line = '%d %s' %(i+1, labels[i])
            print >> outfile, line
            
    # 2.2) If the labels of the nodes are not given
    else:
        for i in xrange(N):
            line = '%d "%d"' %(i+1, i+1)
            print >> outfile, line

    # 3) SAVE THE LINKS
    if adjmatrix[0,0].dtype in [uint8,uint,int8,int]:
        formatstring = '%d %d %d'
    elif adjmatrix[0,0].dtype in [float32,float,float64]:
        formatstring = '%d %d %f'
        
    # 3.1) Save the ARCS if directed
    if directed:
        print >> outfile, '*Arcs'
        for i in xrange(N):
            neighbours = adjmatrix[i].nonzero()[0]
            for j in neighbours:
                line = formatstring %(i+1,j+1,adjmatrix[i,j])
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
                    line = formatstring %(i+1,j+1,adjmatrix[i,j])
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
    if type(nodelist1) == ndarray:
        nodelist1 = list(nodelist1)
    if nodelist2 == None:
        nodelist2 = nodelist1
    else:
        if type(nodelist2) == ndarray:
            nodelist2 = list(nodelist2)
    
    # 1) CREATE LISTS OF INDICES FOR SLICING
    N1 = len(nodelist1)
    N2 = len(nodelist2)
    xindices = []
    for node in nodelist1:
        xindices += [node]*N2
    yindices = nodelist2*N1

    # 2) EXTRACT THE SUBMATRIX AND FINISH
    return adjmatrix[xindices,yindices].reshape(N1,N2)

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
    
    if Reciprocity(adjmatrix) == 1:
        return adjmatrix
    else:
        return 0.5*(adjmatrix + adjmatrix.T)

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
    N = len(adjmatrix)

    laplacianmatrix = identity(N,dtype=adjmatrix.dtype) * adjmatrix.sum(axis=1)
    laplacianmatrix -= adjmatrix
    
    return laplacianmatrix

#def CleanCycles(cyclelist):
#    """Removes repeated cycles from a list of cycles.
    
#    Parameters
#    ----------
#    cyclelist : list
#        A list containing cycles as list or array-like. The cycles shall
#        not contain repeated nodes, e.g. [2,30,5,2] should be expressed
#        as [2,30,5].
    
#    Returns
#    -------
#    cleanlist : list
#        A list with every cycle represented only once. Cycles start with
#        the node with lowest index.
        
#    See Also
#    --------
#    PathsAllinOne : Computes several graph distance and path measures at once.
#    CleanPaths : Removes repeated opposite paths from undirected graphs.
    
#    Notes
#    -----
#    The function is aimed to be used with the output of PathsAllinOne().
    
#    Examples
#    --------
#    >>> cycles = [[3,1,2],[1,2,3],[2,3,1],[5,2],[2,5]]
#    >>> clean = CleanCycles(cycles)
#    >>> clean
#    [[1, 2, 3], [2, 5]]
#    """
#    cleanlist = []
#    for cyc in cyclelist:
#        idx = argmin(cyc)
#        cyc = cyc[idx:] + cyc[:idx]
#        if cyc not in cleanlist:
#            cleanlist.append(cyc)
        
#    return cleanlist

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
def ArrayCompare(a1, a2):
    """Compares whether two arrays are identical or not.
    
    Returns True if all elements are equal, and False if at least
    one element is different.
    
    Parameters
    ----------
    a1 : ndarray
    a2 : ndarray
    """
    
    # 0) SECURITY CHECK
    assert shape(a1) == shape(a2), 'Arrays are not alligned.'
    
    # 1) COMPARE THE ARRAYS
    test = (a1==a2)
    return test.all()

def HammingDistance(array1, array2, normed=False):
    """Computes the Hamming distance between two arrays of same shape.
    
    The Hamming distance counts the fraction of elements that are different
    in two arrays. It is 0 when all the two arrays are identical and 1 when
    the two arrays have no common elements.

    Note
    ----
    In case of binary data, it can be normed to account for the random
    probability of coincidences in two arrays with N1,0 and N2,0 number
    of zero entries and N1,1 and N2,1 number of ones.
    
    Parameters
    ----------
    array1 : ndarray of any shape and dimension
    array2 : ndarray of same shape and dimension as array1
    normed : Boolean, optional
        If True, the function computes the expected Hamming distance between
        the two vectors. The expected value is based on the expected number
        of coincidences that would be found if both arrays were random
        binary vectors with n1 and n2 number of ones respectively (giving
        n - n1 and n - n2 number of zeros).
    
    Returns
    -------
    Hdist : floating point.
        The Hamming distance between the two arrays.
    expHdist : floating point. Optional, only if 'expected = True'
        The expected Hamming distance between the vectors
    """
    # 0) PREPARE FOR CALCULATIONS
    # 0.1) Convert the arrays into rank-1 arrays
    array1 = array1.flatten()
    array2 = array2.flatten()
  
    # 0.2) Security check
    assert len(array1) == len(array2), 'Arrays are not alligned'

    # 1) COUNT THE NUMBER OF COINCIDENCES
    similarity = (array1 == array2)
    Nequal = similarity.sum()

    # 2) COMPUTE THE HAMMING DISTANCE
    Hdist = 1. - float(Nequal)/len(array1)
    
    # 3) RETURN RESULT ACCORDING TO OPTIONS
    # Standard Hamming distance
    if not normed:
        return Hdist
    
    # Normalized Hamming distance
    else:
        # Count the number of ones in the two arrays
        N1 = len(array1.nonzero()[0])
        N2 = len(array2.nonzero()[0])
        
        # Estimate the expected number of random coincidences
        N = len(array1)
        expNc = 1./N * (N1*N2 + (N-N1) * (N-N2))
        
        # The expected Hamming distance
        expHdist = 1. - expNc / N
        
        return Hdist, expHdist

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
        data = array(data)
    
    # 2) Find the minimal non-zero value and return.
    idx = where(data)
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
    ydata, xdata = histogram(data, nbins, range, normed)

    # 1.1) Compute the cumulative sum of the probability
    ydata = ydata.cumsum()

    # 2) RETURN THE RESULTS
    if normed == True:
        norm = 1./ydata[-1]
        ydata *= norm
        
        return xdata[:-1], ydata
    
    elif normed==False:
        return xdata[:-1], ydata

def Factorial(x):
    """Computes the factorial of an integer number.
    """
    # 0) SECURITY CHECK
    assert isinstance(x, int), "'factorial' function only accepts integers"
    
    # 1) COMPUTE THE FACTORIAL
    if x == 0:
        return 1
    
    else:
        fact = x
        for i in xrange(x-1,0,-1):
            fact *= i
            
        return fact

def BinomialCoefficient(n,m):
    """Computes the binomial coefficient of n over m.
    """
    ma = max(n-m,m)
    mi = min(n-m,m)
    enum = 1
    for i in xrange(ma+1,n+1):
        enum *= i
        
    return enum/factorial(mi)

def StdDeviation(data):
    """Returns the mean value and the standard deviation of an array of data.
    It is a simple wrapper using the numpy ufuncs a.mean() and a.dev().
    """
    if type(data) == ndarray:
        return data.mean(), data.std()
    else:
        #print 'Converting data...'
        data = array(data)
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
    # 0) ORDER THE DATA BOTH IN ASCENDING ORDER
    Ndata = len(data)

    if type(data) != ndarray: data = array(data)
    
    minmaxdata = data.copy()
    minmaxdata.sort()

    # 1) FIND THE QUARTILES
    # 1.1) First quartile
    id1 = (Ndata+1) * 0.25
    id1 = int(id1) - 1  # -1 only because of Python starting at 0!!
    Q1 = 0.5* (minmaxdata[id1] + minmaxdata[id1+1])

    # 1.2) The median
    id2 = (Ndata+1) * 0.5
    id2 = int(id2) - 1  # -1 only because of Python starting at 0!!
    Q2 = 0.5* (minmaxdata[id2] + minmaxdata[id2+1])

    # 1.3) Third quartile
    id3 = (Ndata+1) * 0.75
    id3 = int(id3) - 1 # -1 only because of Python starting at 0!!
    Q3 = 0.5* (minmaxdata[id3] + minmaxdata[id3+1])

    return Q1, Q2, Q3

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
    permutationlist = []
    for item in itertools.permutations(data):
        permutationlist.append(item)
    
    return permutationlist

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
    combinationlist = []
    for item in itertools.combinations(data, comblength):
        combinationlist.append(item)
        
    return combinationlist

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
    if float(Ndata)/2 == Ndata/2: isodd = False
    else: isodd = True

    # 1) FIND ALL THE COMBINATIONS UP TO LENGTH len(data)/2
    Nmax = Ndata/2
    setdata = set(data)
    bipartitions = []
    for n in xrange(1,Nmax):
        combinations = AllCombinations(data,n)
        for comb in combinations:
            complementary = tuple(setdata - set(comb))
            bipartitions.append((comb,complementary))

    # 2) ORGANIZSE THE COMBINATIONS INTO BIPARTITIONS
    if isodd:
        combinations = AllCombinations(data,Nmax)
        for comb in combinations:
            complementary = tuple(setdata - set(comb))
            bipartitions.append((comb,complementary))
            
    else:
        combinations = AllCombinations(data,Nmax)
        helperlist = []
        for comb in combinations:
            if set(comb) in helperlist: continue
            # else
            complementary = tuple(setdata - set(comb))
            bipartitions.append((comb,complementary))
            helperlist.append(set(complementary))    
        del helperlist

    return bipartitions

