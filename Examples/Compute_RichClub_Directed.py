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
Here the different cases for DIRECTED graphs are considered.
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
import galib.models
from galib.models import RewireNetwork
from galib.tools import SymmetriseMatrix, LoadFromPajek


##################################################################
# 0) READ THE DATA
currdir = os.getcwd()
dataroot = os.path.join(currdir, 'Data/')
net = loadtxt(dataroot + 'Cat53_cortex.txt', dtype=uint8)
N = len(net)


# 1) COMPUTE THE RICH-CLUB OF THE NETWORKS
# 1.1) Empirical network. Three different cases: input, output & average degrees
# Notice that 'net' is weighted but function RichClub ignores the weights.
inphi = RichClub(net, rctype='indegree')
outphi= RichClub(net, rctype='outdegree')
avphi = RichClub(net, rctype='average')

# 1.2) Rich-club in an ensemble of rewired networks for comparison
nrealiz = 100
prewire = 10
inkmax = len(inphi)
rewphi_in = zeros((nrealiz,inkmax), float)
outkmax = len(outphi)
rewphi_out = zeros((nrealiz,outkmax), float)
avkmax = len(avphi)
rewphi_av = zeros((nrealiz,avkmax), float)
for re in range(nrealiz):
    # Generate a randomly rewired graph (conserving the in-/out-degrees)
    rewnet = RewireNetwork(net, prewire, directed=True)

    # Compute the rich club of the rewired network
    rewphi_in[re] = RichClub(rewnet, rctype='indegree')
    rewphi_out[re] = RichClub(rewnet, rctype='outdegree')
    rewphi_av[re] = RichClub(rewnet, rctype='average')

# Find the average rich-club and deviation for each k value
maxrewphi_in = rewphi_in.max(axis=0)
minrewphi_in = rewphi_in.min(axis=0)
meanrewphi_in = rewphi_in.mean(axis=0)
devrewphi_in = rewphi_in.std(axis=0)

maxrewphi_out = rewphi_out.max(axis=0)
minrewphi_out = rewphi_out.min(axis=0)
meanrewphi_out = rewphi_out.mean(axis=0)
devrewphi_out = rewphi_out.std(axis=0)

maxrewphi_av = rewphi_av.max(axis=0)
minrewphi_av = rewphi_av.min(axis=0)
meanrewphi_av = rewphi_av.mean(axis=0)
devrewphi_av = rewphi_av.std(axis=0)


# 2) PLOT THE RESULTS. REQUIRES MATPLOTLIB!!
colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

# 2.1) Rich-club for the INPUT degrees
plt.figure()
plt.title('Rich-club for input connections', fontsize=14)
# Plot the rich-club of the network
klist = arange(inkmax)
plt.plot(klist, inphi, color=colorlist[1], label='Original network')

# Plot the mean rich-club of rewired networks with errorbars
plt.plot(klist,meanrewphi_in, color=colorlist[0], label='Rewired networks')
plt.plot(klist, maxrewphi_in, '--', color=colorlist[0], zorder=0)
plt.plot(klist, minrewphi_in, '--', color=colorlist[0], zorder=0)

# Eye candy
plt.xlabel('(Input) Node degree', fontsize=14)
plt.ylabel('$\phi$-density', fontsize=14)
plt.legend(loc='upper left', fontsize=10)
plt.grid()


# 2.2) Rich-club for the OUTPUT degrees
plt.figure()
plt.title('Rich-club for output connections', fontsize=14)
# Plot the rich-club of the network
klist = arange(inkmax)
plt.plot(klist, outphi, color=colorlist[1], label='Original network')

# Plot the mean rich-club of rewired networks with errorbars
plt.plot(klist,meanrewphi_out, color=colorlist[0], label='Rewired networks')
plt.plot(klist, maxrewphi_out, '--', color=colorlist[0], zorder=0)
plt.plot(klist, minrewphi_out, '--', color=colorlist[0], zorder=0)

# Eye candy
plt.xlabel('(Output) Node degree', fontsize=14)
plt.ylabel('$\phi$-density', fontsize=14)
plt.legend(loc='upper left', fontsize=10)
plt.grid()


# 2.3) Rich-club for the AVERAGE degrees
plt.figure()
plt.title('Rich-club for average connections', fontsize=14)
# Plot the rich-club of the network
klist = arange(avkmax)
plt.plot(klist, avphi, color=colorlist[1], label='Original network')

# Plot the mean rich-club of rewired networks with errorbars
plt.plot(klist,meanrewphi_av, color=colorlist[0], label='Rewired networks')
plt.plot(klist, maxrewphi_av, '--', color=colorlist[0], zorder=0)
plt.plot(klist, minrewphi_av, '--', color=colorlist[0], zorder=0)

# Eye candy
plt.xlabel('(Averaged) Node degree', fontsize=14)
plt.ylabel('$\phi$-density', fontsize=14)
plt.legend(loc='upper left', fontsize=10)
plt.grid()


plt.show()
