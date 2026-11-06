# -*- coding: utf-8 -*-
# Copyright (c) 2013, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
INPUT / OUTPUT FUNCTIONS
========================

This module contains functions to load and save datasets.

LoadFromPajek
    Reads a network from a text file with Pajek format.
Save2Pajek
    Saves a network into a text file with Pajek format.
LoadLabels
    Reads the labels of nodes from a text file.
SaveLabels
    Saves the labels of nodes into a text file
LoadPartition
    Reads a partition of nodes from a text file.
SavePartition
    Saves a partition of nodes into a text file.

...moduleauthor:: Gorka Zamora-López <gorka@zamora-lopez.xyz>

"""
# Standard library imports
import re
import warnings
# Third party imports
import numpy as np
# Local imports


# DATA I/O FUNCTIONS ###########################################################
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
            warnings.warn( "No labels found to read.", category=RuntimeWarning )

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
                raise ValueError( "Could not find whether network is directed or undirected")
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
        outdtype = np.int64
        try:
            outdtype(aij)
        except ValueError:
            outdtype = np.float64

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
        if adjmatrix[0, 0].dtype in [np.uint8, np.uint, np.int8, int]:
            formatstring = '%d %d %d'
        elif adjmatrix[0, 0].dtype in [np.float16, np.float32, float, np.float64]:
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


#
