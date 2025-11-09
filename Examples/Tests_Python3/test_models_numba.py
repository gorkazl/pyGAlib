# -*- coding: utf-8 -*-
# Copyright (c) 2013, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""In this file I will test all functions in module gamodels_numba.py.
"""

# Standard library imports
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
dataroot = '../Data/'
net, labs = LoadFromPajek(dataroot + 'Cat53_cortex.net', getlabels=True)
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
randnet = RandomGraph_Numba(N,L, directed=False)


# 3) RANDOM DIGRAPHS
dens = 0.1
L = int(dens*N*(N-1))
randnet = RandomGraph_Numba(N,L, directed=True)

time2 = timer()
print( time2 - time1, 'Seconds' )


#
