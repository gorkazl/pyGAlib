"""
In this script we show how to compute the rich-club coefficient of a graph,
including the expected rich-club as quantified out of a set of rewired
networks for comparison. It requires matplotlib only for the visualization.
"""

__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@zamora-lopez.xyz"
__copyright__ = "Copyright 2013-2018"
__license__ = "GPL"
__update__="07/01/2015"

from os.path import join
from numpy import*
import matplotlib.pyplot as plt
#from gatools import*
from galib import RichClub
from gamodels import RewireNetwork


# 0) READ THE DATA
datapath = '/yourpath/GAlib/Data/'
net = loadtxt(join(datapath,'Zachary.txt'), uint8)
N = len(net)

# 1) COMPUTE THE RICH-CLUB OF THE NETWORK
# 1.1) RichClub of the input degrees
phidens = RichClub(net, weightednet=True, rctype='undirected')

# 2) COMPUTE THE RICH-CLUB IN AN ENSEMBLE OF REWIRED NETWORKS FOR COMPARISON
nrealiz = 100
prewire = 10
kmax = len(phidens)
rewphidens = zeros((nrealiz,kmax), float)
for re in xrange(nrealiz):
    # 2.1) Rewire the network
    rewnet = RewireNetwork(net, prewire, directed=False)

    # 2.2) Compute the rich club of the rewired network
    rewphidens[re] = RichClub(rewnet, rctype='undirected')

# 2.3) Find the average rich-club and deviation for each k value
meanrewphi = rewphidens.mean(axis=0)
devrewphi = rewphidens.std(axis=0)


# 3) PLOT THE RESULTS. REQUIRES MATPLOTLIB!!
plt.figure()

# 3.1) Plot the rich-club of the network
xdata = arange(kmax)
plt.plot(xdata, phidens, label='Original network')

# 3.2) Plot the mean rich-club of rewired networks with errorbars
plt.errorbar(xdata,meanrewphi,devrewphi, label='Rewired networks')

# 3.3) Eye candy
plt.xlabel('degree', fontsize=16)
plt.ylabel('$\phi$-density', fontsize=16)
plt.legend(loc='upper left')


plt.show()
