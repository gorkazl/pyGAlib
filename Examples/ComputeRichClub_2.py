"""
In this script we show how to compute the rich-club coefficient of a graph,
including the expected rich-club as quantified out of a set of rewired 
networks for comparison. Here we illustrate the use of RichClub() function 
in case of directed networks, for which different manners to compute the
rich-club are possible. It requires matplotlib only for the visualization.
"""

__author__ = "Gorka Zamora-Lopez" 
__email__ = "Gorka.zamora@ymail.com"
__copyright__ = "Copyright 2013-2015"
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
net = loadtxt(join(datapath,'Cat53_cortex.txt'), uint8)
N = len(net)

# 1) COMPUTE THE RICH-CLUB OF THE NETWORK
# 1.1) RichClub of the input degrees
inphi = RichClub(net, weightednet=True, rctype='indegree')
outphi= RichClub(net, weightednet=True, rctype='outdegree')
avphi = RichClub(net, weightednet=True, rctype='average')

# 2) COMPUTE THE RICH-CLUB IN AN ENSEMBLE OF REWIRED NETWORKS FOR COMPARISON
nrealiz = 5
prewire = 10
inkmax = len(inphi)
inrewphi = zeros((nrealiz,inkmax), float)
outkmax = len(outphi)
outrewphi = zeros((nrealiz,outkmax), float)
avkmax = len(avphi)
avrewphi = zeros((nrealiz,avkmax), float)
for re in xrange(nrealiz):
    # 2.1) Rewire the network
    rewnet = RewireNetwork(net, prewire, directed=True)
    
    # 2.2) Compute the rich club of the rewired network
    inrewphi[re] = RichClub(rewnet, rctype='indegree')
    outrewphi[re] = RichClub(rewnet, rctype='outdegree')
    avrewphi[re] = RichClub(rewnet, rctype='average')

# 2.3) Find the average rich-club and deviation for each k value
meaninrewphi = inrewphi.mean(axis=0)
devinrewphi = inrewphi.std(axis=0)

meanoutrewphi = outrewphi.mean(axis=0)
devoutrewphi = outrewphi.std(axis=0)

meanavrewphi = avrewphi.mean(axis=0)
devavrewphi = avrewphi.std(axis=0)


# 3) PLOT THE RESULTS. REQUIRES MATPLOTLIB!!
## RICH-CLUB FOR THE INPUT DEGREES
plt.figure()
kmax = max(inkmax,outkmax,avkmax)

# 3.1) Plot the rich-club of the network
xdata = arange(inkmax)
plt.plot(xdata, inphi, label='Original network')

# 3.2) Plot the mean rich-club of rewired networks with errorbars
plt.errorbar(xdata,meaninrewphi,devinrewphi, label='Rewired networks')

# 3.3) Set add-ons of the figure
plt.xlim(0,kmax)
plt.ylim(0,1)

plt.title('Rich-club of input degrees')
plt.xlabel('degree', fontsize=16)
plt.ylabel('$\phi$-density', fontsize=16)
plt.legend(loc='upper left')


## RICH-CLUB FOR THE OUTPUT DEGREES
plt.figure()

# 3.1) Plot the rich-club of the network
xdata = arange(outkmax)
plt.plot(xdata, outphi, label='Original network')

# 3.2) Plot the mean rich-club of rewired networks with errorbars
plt.errorbar(xdata,meanoutrewphi,devoutrewphi, label='Rewired networks')

# 3.3) Set add-ons of the figure
plt.xlim(0,kmax)
plt.ylim(0,1)

plt.title('Rich-club of output degrees')
plt.xlabel('degree', fontsize=16)
plt.ylabel('$\phi$-density', fontsize=16)
plt.legend(loc='upper left')


## RICH-CLUB FOR THE AVERAGE DEGREES
plt.figure()

# 3.1) Plot the rich-club of the network
xdata = arange(avkmax)
plt.plot(xdata, avphi, label='Original network')

# 3.2) Plot the mean rich-club of rewired networks with errorbars
plt.errorbar(xdata,meanavrewphi,devavrewphi, label='Rewired networks')

# 3.3) Set add-ons of the figure
plt.xlim(0,kmax)
plt.ylim(0,1)

plt.title('Rich-club of average degrees')
plt.xlabel('degree', fontsize=16)
plt.ylabel('$\phi$-density', fontsize=16)
plt.legend(loc='upper left')


plt.show()


