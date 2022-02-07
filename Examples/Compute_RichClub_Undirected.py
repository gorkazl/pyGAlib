# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2022, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
This script computes the rich-club coefficient of a graph. It also does so for
an ensemble of rewired networks conserving the degree distribution.
UNDIRECTED graphs considered only.
"""
from __future__ import division, print_function

__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@zamora-lopez.xyz"
__copyright__ = "Copyright 2013-2022"
__license__ = "Apache License 2.0"
__update__="07/02/2022"

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
from galib.tools import SymmetriseMatrix, LoadFromPajek


##################################################################
# 0) READ THE DATA
currdir = os.getcwd()
dataroot = os.path.join(currdir, 'Data/')
net = loadtxt(dataroot + 'Zachary.txt', dtype=uint8)
net = LoadFromPajek(dataroot + 'Dolphins.net', getlabels=False)
net, labels = LoadFromPajek(dataroot + 'LesMiserables.net', getlabels=True)
# net = loadtxt(dataroot + 'Cat53_cortex.txt', dtype=uint8)
# net = SymmetriseMatrix(net)
N = len(net)


# 1) COMPUTE THE RICH-CLUB OF THE NETWORKS
# 1.1) Rich-club of the empirical network.
# Notice that 'net' is weighted but function RichClub ignores the weights.
rcdens = RichClub(net, rctype='undirected')

# 1.2) Rich-club in an ensemble of rewired networks for comparison
nrealiz = 100
prewire = 10
kmax = len(rcdens)
rewrcdens = zeros((nrealiz,kmax), float)
for re in range(nrealiz):
    # Generate a randomly rewired graph (conserving the degrees)
    rewnet = RewireNetwork(net, prewire, directed=False)

    # Compute the rich club of the rewired graph
    rewrcdens[re] = RichClub(rewnet, rctype='undirected')

# Find the average rich-club and deviation for each k value
maxrewphi = rewrcdens.max(axis=0)
minrewphi = rewrcdens.min(axis=0)
meanrewphi = rewrcdens.mean(axis=0)
devrewphi = rewrcdens.std(axis=0)


# 2) PLOT THE RESULTS. REQUIRES MATPLOTLIB!!
plt.figure()
colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

# 2.1) Plot the rich-club of the network
klist = arange(kmax)
plt.plot(klist, rcdens, color=colorlist[1], label='Original network')

# 2.2) Plot the mean rich-club of rewired networks with errorbars
# plt.errorbar(klist,meanrewphi,devrewphi, color=colorlist[1], label='Rewired networks')
plt.plot(klist,meanrewphi, color=colorlist[0], label='Rewired networks')
plt.plot(klist, maxrewphi, '--', color=colorlist[0], zorder=0)
plt.plot(klist, minrewphi, '--', color=colorlist[0], zorder=0)

# 2.3) Eye candy
plt.xlabel('Node degree', fontsize=14)
plt.ylabel('$\phi$-density', fontsize=14)
plt.legend(loc='upper left', fontsize=10)
plt.grid()


plt.show()
