# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
In this script we show how to convert datasets
"""
from __future__ import division, print_function

__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@zamora-lopez.xyz"
__copyright__ = "Copyright 2013-2019"
__license__ = "Apache Lincense 2.0"
__update__="13/06/2019"

# Standard library imports
import os, os.path
# Third party imports
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from numpy import*
# Personal libraries
from galib import RichClub
from galib.models import RewireNetwork
from galib.tools import LoadFromPajek, Save2Pajek, LoadLabels, SaveLabels


##################################################################
# 1) READ SOME DATA IN PAJEK FORMAT AND SAVE THE ADJACENCY MATRIX
currdir = os.getcwd()
dataroot = os.path.join(currdir, 'Data/')

# 1.1) Read the data splitting the adjacency matrix and the labels
fname = 'LesMiserables.net'
net, labels = LoadFromPajek(dataroot + fname, getlabels=True)

# 1.2) Save the adjacency matrix as a text file
# Give '%d' formatter for integer data, '%f' for real valued data
outfname1 = 'spam_LesMiserables.txt'
savetxt(dataroot + outfname1, net, fmt='%d')

# 1.3) Save the adjacency matrix in a numpy binary file
outfname2 = 'spam_LesMiserables.npy'
save(dataroot + outfname2, net)

# 1.4) Save the labels in an independent text file
outfname3 = 'spam_LesMiserables_labels.txt'
SaveLabels(dataroot + outfname3, labels)


# 2) READ AN ADJACENCY MATRIX AND SAVE IT AS PAJEK FORMAT
# 2.1) Read the adjacency matrix
fname = 'Cat53_cortex.txt'
net = loadtxt(dataroot + fname, dtype=uint8)

# 2.2) Read the labels from a file
fname2 = 'Cat53_labels.txt'
netlabels = LoadLabels(dataroot + fname2)

# 2.3) Save the network and labels into a Pajek formatted text file
outfname = 'spam_Cat53_cortex.net'
Save2Pajek(dataroot + outfname, net, labels=netlabels, directed=True)
