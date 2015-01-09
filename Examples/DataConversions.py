"""
In this script we show how to convert datasets
"""

__author__ = "Gorka Zamora-Lopez" 
__email__ = "Gorka.zamora@ymail.com"
__copyright__ = "Copyright 2013-2015"
__license__ = "GPL"
__update__="07/01/2015"

from os.path import join
from numpy import*
from gatools import LoadFromPajek, Save2Pajek, LoadLabels, SaveLabels


# 1) READ SOME DATA IN PAJEK FORMAT AND SAVE THE ADJACENCY MATRIX
datapath = 'Data/'

# 1.1) Read the data splitting the adjacency matrix and the labels
fname = 'LesMiserables.net'
net, labels = LoadFromPajek(join(datapath,fname), getlabels=True)

# 1.2) Save the adjacency matrix as a text file
# Give '%d' formatter for integer data, '%f' for real valued data
outfname1 = 'LesMiserables.txt'
savetxt(join(datapath,outfname1), net, fmt='%d')

# 1.3) Save the adjacency matrix in a numpy binary file
outfname2 = 'LesMiserables.npy'
save(join(datapath,outfname2), net)

# 1.4) Save the labels in an independent text file
outfname3 = 'LesMiserables_labels.txt'
SaveLabels(join(datapath,outfname3), labels)


# 2) READ AN ADJACENCY MATRIX AND SAVE IT AS PAJEK FORMAT
datapath = 'Data/'

# 2.1) Read the adjacency matrix
fname = 'Cat53_cortex.txt'
net = loadtxt(join(datapath,fname), dtype=uint8)

# 2.2) Read the labels from a file
fname2 = 'Cat53_labels.txt'
netlabels = LoadLabels(join(datapath,fname2))

# 2.3) Save the network and labels into a Pajek formatted text file
outfname = 'Cat53_cortex.net'
Save2Pajek(join(datapath,outfname), net, labels=netlabels, directed=True)



