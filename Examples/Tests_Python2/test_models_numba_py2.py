# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2018, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""In this file I will test all functions in the Python3 version of GAlib
(the module gamodels_numba.py) and make sure they work.
"""

# Standard library imports
import os, os.path
from timeit import default_timer as timer
# Third party imports
import matplotlib.pyplot as plt
from numpy import*
# Personal libraries
from galib.tools import LoadFromPajek
from galib.models import RandomGraph
from galib.models_numba import RandomGraph_Numba

##################################################################
# 0) READ THE DATA
currdir = os.getcwd()
currdir = os.path.split(currdir)[0]
dataroot = os.path.join(currdir, 'Data/')
net, labs = LoadFromPajek(dataroot + 'Cat53_cortex.net', True)
netsym = 0.5*(net+net.T)
N = len(net)

# Define the partition
visual = arange(16)
audit = arange(16,23)
somatomotor = arange(23,39)
frontolimbic = arange(39, 53)
partition = [visual,audit,somatomotor,frontolimbic]
ncoms = len(partition)

time1 = timer()
N = 5000
# # 1) RING MODELS
# latt = Lattice1D(N,5)
# latt2 = Lattice1D_FixLinks(N,550)
# wsnet = WattsStrogatzGraph(N,5, 0.05)


# 2) RANDOM GRAPHS
dens = 0.1
L = int(0.5*dens*N*(N-1))
# ernet = ErdosRenyiGraph(N, dens)
# randnet = RandomGraph(N,L)
randnet = RandomGraph_Numba(N,L)
# sfnet = ScaleFreeGraph(N,dens, 3)


# 3) RANDOM DIGRAPHS
dens = 0.1
L = int(dens*N*(N-1))
# ernet = ErdosRenyiGraph(N, dens, directed=True)
# randnet = RandomGraph(N,L, directed=True)
randnet = RandomGraph_Numba(N,L, directed=True)
# banet = BarabasiAlbertGraph(N,5)
# sfnet = ScaleFreeGraph(N,dens, 3, directed=True)

#
# # 4) REWIRING ALGORITHMS
# # rewnet = RewireNetwork(net)
# # rewnet2 = ModularityPreservingGraph(net, partition)
# rewnet = RewireNetwork(netsym)
# rewnet2 = ModularityPreservingGraph(netsym, partition)
# Qorig = Modularity(netsym.astype(bool), partition)
# Qrew = Modularity(rewnet2.astype(bool), partition)
# print 'Modularity', Qorig, Qrew
#
#
# # 5) MODULAR & HIERARCHICAL NETS
# modnet, newpart = ModularHeterogeneousGraph([100,200,400], [0.2,0.1,0.05], 0.02, directed=True)
# # testpart = HMpartition([100,200,400])
# hmrnet = HMRandomGraph([2,3,100], [2,10,20], directed=True)
# hmcnet = HMCentralizedGraph([2,3,100], [2,10,20], [1.8,2.5,3.0], directed=False)


time2 = timer()
print time2 - time1, 'Seconds'

#
