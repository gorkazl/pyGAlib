# -*- coding: utf-8 -*-
# Copyright (c) 2013, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
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

# Standard library imports
# Third party imports
import matplotlib.pyplot as plt
import numpy as np
# Local libraries
import galib


##################################################################
# 0) READ THE DATA
dataroot = 'Data/'
# net = np.loadtxt(dataroot + 'Zachary.txt', dtype=np.uint8)
# net, labels = galib.tools.LoadFromPajek(dataroot + 'LesMiserables.net', getlabels=True)
net = np.loadtxt(dataroot + 'Cat53_cortex_sym.txt').astype(np.uint8)

# Get some basic properties of the network
N = len(net)
L = int( 0.5 * net.astype(bool).sum() )
deg = galib.Degree(net)


# 1) COMPUTE THE k-DENSITY OF THE NETWORKS
# 1.1) k-density of the empirical graph
# Notice that 'net' can be weighted but function RichClub ignores the weights.
kdens = galib.k_Density(net, rctype='undirected')

# 1.2) Surrogates - Compute k-density in RANDOM graphs
nrealiz = 100

kmax = len(kdens)
klist = np.arange(kmax)

deg_rand = np.zeros((nrealiz,N), np.uint64)
kdens_rand = np.zeros((nrealiz,kmax), np.float64)
for re in range(nrealiz):
    # Generate the rewired surrogate graph (conserving the degrees)
    randnet = galib.RandomGraph(N,L, directed=False)
    deg_rand[re] = galib.Degree(randnet)

    # Compute its k-density
    spam = galib.k_Density(randnet, rctype='undirected')
    kdens_rand[re,:len(spam)] = spam

# Find the average rich-club and deviation for each k value
maxphi_rand = kdens_rand.max(axis=0)
minphi_rand = kdens_rand.min(axis=0)
meanphi_rand = kdens_rand.mean(axis=0)


# 1.3) Surrogates - Compute k-density in REWIRED graphs
nrealiz = 100
prewire = 10
deg_rew = np.zeros((nrealiz,N), np.uint64)
kdens_rew = np.zeros((nrealiz,kmax), np.float64)
for re in range(nrealiz):
    # Generate the rewired surrogate graph (conserving the degrees)
    rewnet = galib.RewireNetwork(net, prewire, directed=False)
    deg_rew[re] = galib.Degree(rewnet)

    # Compute its k-density
    kdens_rew[re] = galib.k_Density(rewnet, rctype='undirected')

# Find the average rich-club and deviation for each k value
maxphi_rew = kdens_rew.max(axis=0)
minphi_rew = kdens_rew.min(axis=0)
meanphi_rew = kdens_rew.mean(axis=0)



# 2) PLOT THE RESULTS. REQUIRES MATPLOTLIB!!
plt.figure(figsize=(6.4,6.4))

# 2.1) Plot the degree distributions
plt.subplot(2,1,1)
plt.hist(deg,              bins=N, range=[0,N], width=0.8, color='tab:red', alpha=0.8, align='mid', density=True , label='Empirical', zorder=0)
plt.hist(deg_rand.ravel(), bins=N, range=[0,N], width=0.6, color='gray', alpha=0.5, align='mid', density=True, label='Random Graphs', zorder=2 )
# plt.hist(deg_rew.ravel(),  bins=N, range=[0,N], width=0.6, color='tab:blue', alpha=0.8, align='mid', density=True, label='Rewired Grahs', zorder=1 )
plt.xlim(-1,kmax+1)
plt.ylabel('Probability')
plt.legend(loc='upper right')


# 2.2) Plot the k-densities
plt.subplot(2,1,2)
plt.plot(klist, kdens, '.-', color='tab:red', label='Empirical')
plt.plot(klist, meanphi_rand, color='gray', label='Random graphs')
plt.fill_between(klist, minphi_rand, maxphi_rand, color='gray', alpha=0.3, linewidth=1)
plt.plot(klist, meanphi_rew, color='tab:blue', label='Rewired graphs')
plt.fill_between(klist, minphi_rew, maxphi_rew, color='tab:blue', alpha=0.3)
plt.xlim(-1,kmax+1)
plt.ylim(-0.02,1.05)
plt.xlabel('Node Degree')
plt.ylabel('k-density')
plt.grid(axis='y', ls='-.', lw=0.5)
plt.legend(loc='upper left')

plt.tight_layout()



# 3) IDENTIFY THE SET OF RICH-CLUB HUBS, WITH k-DENSITY LARGER THAN A THRESHOLD
threshold = 0.8
# Identify the degree at which k-density overcomes the given threshold

if kdens.max() < threshold:
    print( f"No rich-club found for k-density > {threshold:1.3f}" )
else:
    k_crit = np.where(kdens >= threshold)[0][0]

    # Identify indices of the nodes forming the rich-club
    richclub = np.where(deg >= k_crit)[0]
    print( f'List of Rich-Club nodes with k-density >= {threshold:1.3f}' )
    print( richclub )


plt.show()
##
