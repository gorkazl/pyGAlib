# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2022, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
This script is an example of how to compute and plot the roles of nodes
according to the partition of a network into communities.

Citation
--------

F. Klimm, J. Borge-Holthoefer, N. Wessel, J. Kurths, G. Zamora-Lopez (2014)
Individual node's contribution to the mesoscale of complex networks, New J.
Phys. 16, 125006.

"""
from __future__ import division, print_function

__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@zamora-lopez.xyz"
__copyright__ = "Copyright 2013-2022"
__license__ = "Apache Lincese 2.0"
__update__="07/02/2022"

# Standard library imports
import os, os.path
# Third party imports
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from numpy import*
# Personal libraries
from galib import*
from galib.tools import SymmetriseMatrix

def MaxDispersionCurve(d, M):
    """This function is defined only to be plotted in the dispersion vs
    participation plot.
    """
    M = float(M)
    x = 0.5 * (2.0 - d)

    return 1.0 - sqrt( M/(M-1.0) * ( x**2 + (1.0 - x)**2 - 1.0/M ) )


############################################################################
# 0) READ THE DATA OF THE CORTICAL NETWORK OF THE CAT
currdir = os.getcwd()
dataroot = os.path.join(currdir, 'Data/')
net = loadtxt(dataroot + 'Cat53_cortex.txt', uint8)
net = SymmetriseMatrix(net)
N = len(net)

# Define the partition of cat cortical network into functional subdivision
visual = arange(16)
auditory = arange(16,23)
somatomotor = arange(23,39)
frontolimbic = arange(39, 53)
partition = [visual, auditory, somatomotor, frontolimbic]
ncoms = len(partition)


# 1) COMPUTE THE FOUR PARAMETERS DEFINIG THE ROLE OF EACH NODE --
# global hubness, local hubness, node participation and node dispersion
ghubness, lhubness, pindex, dindex = RolesNodes(net, partition)


# 2) PLOT THE RESULTS
# 2.0) Some helper data for the plots
# Define colors of the nodes
colorlist = ['black', 'red', 'green', 'blue']
colornodes = []
for n in range(ncoms):
    com = partition[n]
    colornodes += len(com)*[colorlist[n]]

plt.figure(figsize=(10,3.5))

# 2.1) Plot local hubness vs global hubness
plt.subplot(1,3,1)

# Plot the data
plt.scatter(lhubness,ghubness, color=colornodes, s=20)

# Draw the axis (look for other manner!)
plt.plot((-100,100),(0,0), color='gray', zorder=0)
plt.plot((0,0),(-100,100), color='gray', zorder=0)

# Draw the significance lines
plt.plot((-100,100),(2.5,2.5), '--', color='gray', zorder=0)
plt.plot((-100,100),(-2.5,-2.5), '--', color='gray', zorder=0)
plt.plot((2.5,2.5),(-100,100), '--', color='gray', zorder=0)
plt.plot((-2.5,-2.5),(-100,100), '--', color='gray', zorder=0)

# Axis properties
#plt.xlim(floor(lhubness.min()), floor(lhubness.max()+1))
#plt.ylim(floor(ghubness.min()), floor(ghubness.max()+1))
plt.xlim(-5,5)
plt.ylim(-6,6)
plt.xlabel('Local hubness', fontsize=14)
plt.ylabel('Global hubness', fontsize=14)


# 2.2) Plot node participation vs. global hubness
plt.subplot(1,3,2)

# Plot the data
plt.scatter(pindex,ghubness, color=colornodes, s=20)

# Draw the axis (look for other manner!)
plt.plot((-100,100),(0,0), color='gray', zorder=0)
#plt.plot((0,0),(-100,100), color='gray', zorder=0)

# Draw the significance lines
plt.plot((-100,100),(2.5,2.5), '--', color='gray', zorder=0)
plt.plot((-100,100),(-2.5,-2.5), '--', color='gray', zorder=0)
plt.plot((0.3333,0.3333),(-100,100), '--', color='gray', zorder=0)
plt.plot((0.6666,0.6666),(-100,100), '--', color='gray', zorder=0)

# Axis properties
plt.xlim(-0.01,1.01)
plt.ylim(-6,6)
plt.xlabel('Participation Index', fontsize=14)
plt.ylabel('Global hubness', fontsize=14)

# 2.3) Plot node dispersion vs. node participation
# Compute the background limit curves for the dindex vs. pindex plot
ddata = arange(-0.1,1.001,0.001)
pdata = MaxDispersionCurve(ddata, 4)

plt.subplot(1,3,3)

# Plot the background lines
plt.plot(pdata, ddata, color='gray', zorder=0)
plt.plot((0,1),(0,1), color='gray', zorder=0)

# Include the shades
plt.fill_between((0,1.2), (-0.5,-0.5), (0,1.2), color='gray', alpha=0.2)
plt.fill_between(pdata, ddata, 1.1*ones(len(ddata)), color='gray', alpha=0.2)

# Plot the data
plt.scatter(pindex,dindex, color=colornodes, s=20)
plt.xlabel('Participation Index', fontsize=14)
plt.ylabel('Dispersion Index', fontsize=14)
plt.xlim(-0.01,1.01)
plt.ylim(-0.01,1.01)
plt.grid()

plt.tight_layout()
plt.show()
