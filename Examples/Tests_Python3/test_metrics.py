# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2022, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""In this file I will test all functions in the Python3 version fo GAlib
library and make sure they work.
This file also tests for gametrics_numba.py module. So far, I only have two
functions there so... did not make much sense to have a separate test file.
"""

# Standard library imports
import os, os.path
from timeit import default_timer as timer
# Third party imports
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from numpy import*
# Personal libraries
from galib import*
from galib.tools import*
from galib.metrics_numba import*


##################################################################
# 0) READ THE DATA
currdir = os.getcwd()
currdir = os.path.split(currdir)[0]
dataroot = os.path.join(currdir, 'Data/')
net, labs = LoadFromPajek(dataroot + 'Cat53_cortex.net', True)
N = len(net)

time1 = timer()
# 1) USE galib TO ANALYSE THE NETWORK
dens = Density(net)

ink, outk = Degree(net, True)
ins, outs = Intensity(net, True)

recip = Reciprocity(net)
print('Reciprocity', recip)

kr, kminus, kplus = ReciprocalDegree(net, normed=True)

avkneigh = AvNeighboursDegree(net, knntype='inin')

netsym = (net+net.T).astype('bool')
C, Cnodes = Clustering(netsym)
print('Clustering', C)

kdens = RichClub(netsym)

mimat = MatchingIndex(net)
mimat_numba = MatchingIndex_Numba(net)

dij = FloydWarshall(net)
dij_numba = FloydWarshall_Numba(net)
avpathlen = ( dij.sum() - dij.trace() ) / (N*(N-1))
print('Average pathlength:', avpathlen)

dij, bcnodes, paths, cycles = PathsAllinOne(net)
print('Average pathlength:', avpathlen)

allshortest = ShortestPaths(net, 0,1, dij[0,1])


# DEFINE THE PARTITION
visual = arange(16)
audit = arange(16,23)
somatomotor = arange(23,39)
frontolimbic = arange(39, 53)
partition = [visual,audit,somatomotor,frontolimbic]
ncoms = len(partition)

assmat = AssortativityMatrix(net, partition)
Qnet = Modularity(net.astype(bool), partition)
Qsym = Modularity(netsym, partition)

print('Modularity (dir)', Qnet)
print('Modularity (undir)', Qsym)

kcore9 = K_Core(netsym, kmin=9)
kshells = K_Shells(netsym)

comnode = zeros(N,int)
for c in range(ncoms):
    for node in partition[c]:
        comnode[node] = c


# 1) CALCULATE THE COMMUNITY PARTICIPATION VECTOR FOR EACH NODE
pmatrix = ParticipationMatrix(netsym,partition)


# 2) ROLES OF GUIMERA AND AMARAL
# 2.1) The local hubness of the nodes
zlocal = Hubness_GA(pmatrix, partition)

# 2.2) The participation index
pindex_ga = ParticipationIndex_GA(pmatrix)


# 3) NORMALIZED ROLES
# 3.1) Normalized hubness
degree = Degree(net)
normdeg = degree.astype(float) / N

# 3.2) The participation index
pindex = NodeParticipation(netsym,partition)

ghubness, lhubness, pindex, dindex = RolesNodes(netsym, partition)

time2 = timer()
print( time2 - time1, 'Seconds' )


# 4) PLOT SOME RESULTS
# 4.1) The old indices
plt.figure()
plt.plot(pindex_ga, zlocal, '.')
plt.xlim(0,1)
plt.grid()

# 4.2) The new indices
plt.figure()
plt.plot(pindex, ghubness, '.')
plt.xlim(0,1)
plt.grid()

plt.figure()
plt.plot(lhubness, ghubness, '.')
plt.grid()

plt.figure()
plt.plot(pindex, dindex, '.')
plt.grid()

plt.show()




#
