"""
In this script we show an example of how to compute the roles of nodes
according to the partition of a network into communities.
"""

from os.path import join
from numpy import*
import matplotlib.pyplot as plt
from gatools import SymmetriseMatrix
from galib import *
from gamodels import *

# 0) CREATE A RANDOM MODULAR / HIERARCHICAL NETWORK
# We create a network composed of 3 hierarchical levels. The network is
# divided in 4 communites, each of them has for communities of 16 nodes.
# In total the network has 256 nodes.

net = HMRandomNetwork([4,4,16],[1,3,10])
#net = HMRandomNetwork([4,4,16],[10,3,10])
#net = HMRandomNetwork([4,4,16],[10,10,10])

# For comparison, try with a random graph by uncomenting the following line
#net = ErdosRenyiGraph(256,0.1)
N = len(net)

# Define the partition for the cat cortical network
partition = [arange(64), arange(64,128), arange(128,192), arange(192,256)]
ncoms = len(partition)


# 1) CALCULATE THE PARTICIPATION MATRIX
# pmatrix[i,s] is the number of neighbours of node i in community s
pmatrix = ParticipationMatrix(net, partition)


# 2) COMPUTE THE ROLES OF THE NODES
# 2.1) Normalized global hubness of the nodes
degree = Degree(net)
normdeg = degree.astype(float) / N

# 2.2) The participation index of the nodes
for re in xrange(1000):
    pindex = ParticipationIndex(pmatrix,partition)


# 3) FOR COMPARISON, COMPUTE THE ROLES ORIGINALLY GIVEN BY GUIMERA AND AMARAL
# 3.1) The local hubness of the nodes
zlocal = LocalHubness_GA(pmatrix, partition)

# 3.2) The participation index of the nodes
pindex_GA = ParticipationIndex_GA(pmatrix)


# 4) PLOT THE ROLES IN BOTH FRAMEWORKS
# NOTE!! Plotting requires to have Matplotlib installed.
# Otherwise, comment the rest of the script

# Colors for the nodes according to their community
colornodes = ['black']*64 + ['red']*64 + ['green']*64 + ['blue']*64

# 4.1) Plot the roles of the nodes
plt.figure()

plt.scatter(pindex, normdeg, s=30, color=colornodes)
plt.xlabel('Participation index, $P_i$', fontsize=18)
plt.xlim(0,1)
plt.ylabel('Normalized degree, $k_i / N$', fontsize=18)
plt.ylim(0,1)

plt.grid()

# 4.2) Plot the roles in the z-Pi plane as done by Guimera & Amaral
plt.figure()

plt.scatter(pindex_GA, zlocal, s=30, color=colornodes)
plt.xlabel('Participation index of GA, $P_i$', fontsize=18)
plt.xlim(0,1)
plt.ylabel('Within-module degree, $Z$', fontsize=18)
plt.ylim(-3,6)

# Plot the delimiter lines
plt.plot((0,1),(2.5,2.5), 'k')
plt.plot((0.1,0.1), (-10,2.5), 'k')
plt.plot((0.65,0.65), (-10,2.5), 'k')
plt.plot((0.8,0.8), (-10,2.5), 'k')
plt.plot((0.3,0.3), (2.5, 10), 'k')
plt.plot((0.75,0.75), (2.5,10), 'k')


plt.show()

