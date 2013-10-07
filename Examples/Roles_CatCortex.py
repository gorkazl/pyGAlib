"""
In this script we show an example of how to compute the roles of nodes
according to the partition of a network into communities.
"""

from os.path import join
from numpy import*
import matplotlib.pyplot as plt
from gatools import SymmetriseMatrix
from galib import *

# 0) READ THE DATA OF THE CORTICAL NETWORK OF THE CAT
datapath = '/yourpath/GAlib/Data/'
net = loadtxt(join(datapath, 'Cat53_cortex.txt'), uint8)
net = SymmetriseMatrix(net)
N = len(net)

# Define the partition for the cat cortical network
visual = arange(16)
auditory = arange(16,23)
somatomotor = arange(23,39)
frontolimbic = arange(39, 53)
partition = [visual, auditory, somatomotor, frontolimbic]
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
colornodes = ['black']*16 + ['red']*7 + ['green']*16 + ['blue']*14

# 4.1) Plot the roles of the nodes
plt.figure()

plt.scatter(pindex, normdeg, s=50, color=colornodes)
plt.xlabel('Participation index, $P_i$', fontsize=18)
plt.xlim(0,1)
plt.ylabel('Normalized degree, $k_i / N$', fontsize=18)
plt.ylim(0,1)

plt.grid()

# 4.2) Plot the roles in the z-Pi plane as done by Guimera & Amaral
plt.figure()

plt.scatter(pindex_GA, zlocal, s=50, color=colornodes)
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

