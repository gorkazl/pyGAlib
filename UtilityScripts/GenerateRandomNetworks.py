"""
This script is to easily generate ensembles of random networks even for
users with no or little programming experience. After setting of some 
parameters to personalize the results, the script generates the networks
and saves them into files for later use. The main choices are:

NETWORK PARAMETERS
- N : Size of the network (number of nodes).
- realiz : Number of networks to be generated.
- directed : True to generate directed network, False for undirected. 
- model : 'ER' or 'RG', type of random networks desired.
The Erdos-Renyi (ER) model generates random networks with given link
probability, p. Therefore individual realizations will have slightly
different number of links.
The 'RG' option generates random networks all with the same number of links.
- p : Link probability for the 'ER' case.
- L : Number of links in the 'RG' case.

OUTPUT PARAMETERS
- outpath : Path where the resulting networks will be saved.
- filename : Basename for the output files. Final name will be
basename_reX.dat where X is the realization number.
- filetype : Choice to save networks. 'text' saves the adjacency matrices
as text files with a .dat extension. 'binary' saves the adjacency matrices
in the numpy binary format with a .npy extension. 'pajek' saves the networks
into text files with the Pajek format and extension .net.

RUNNING THE SCRIPT
After personalizing the network and output options in the script, open a
terminal window and move to the folder containing the script, e.g. by typing:

> cd ~/GraphAnalysisLibrary/HelperScripts

then type the following command and wait

> python GenerateRandomNetworks.py

"""

from os.path import join
from numpy import save, savetxt, arange
from gatools import Save2Pajek
from gamodels import ErdosRenyiGraph, RandomGraph


### PARAMETERS TO BE MODIFIED BY THE USER ##################################
### NETWORK OPTIONS
# Size of the network, number of nodes
N = 1000

# Number of realizations
realiz = 100

# Directedness. 'True' to generate directed networks, 'False' for undirected
directed = False

# Random network model.
#'ER' for Erdos-Renyi model, or
#'RG' for random networks with given number of links.
netmodel = 'ER'
if netmodel == 'ER':
    # Link probability of the Erdos-Renyi model
    p = 0.1
elif netmodel == 'RG':
    # Number of links in the random graph model
    L = int(0.1*N*(N-1))


#### OUTPUT OPTIONS
# Path where networks will be saved. The folder must already exist.
outpath = '/yourpath/temp/'

# Base name for the files. Final name will be 'basename_reX', where X is
# the realization number
filename = 'RandomNetwork'

# Type of file to be saved.
# - 'text' saves adjacency matrices as text files with .dat extension.
# - 'binary' saves adjacency matrices in numpy binary format with .npy extension.
# - 'pajek' saves networks into Pajek-readable text files with extension .net.
filetype = 'binary'



############################################################################
### UNLESS CONFIDENT, PLEASE DO NOT MODIFY THE CODE FROM THIS LINE.
# 0) GIVE SOME FEEDBACK TO THE USER
print '- Number of nodes:', N

if netmodel == 'ER':
    print '- Erdos-Renyi networks with p =', p
elif netmodel == 'RG':
    if directed:
        dens = float(L) / (N*(N-1))
    else:
        dens = 2*float(L) / (N*(N-1))
    print '- Random networks with', L, 'links (density = %0.5f)' %dens

if directed: print '- Directed'
else: print '- Undirected'


# 1) DO THE WORK
print '\nGenerating and saving %d random networks...' %realiz

for re in xrange(int(realiz)):
    # 1.0) Feedback the state of the calculation
    if realiz <= 200:
        if re in arange(25,realiz,25):
            print re
    else:
        if re in arange(100,realiz,100):
            print re

    # 1.1) Create the network
    if netmodel == 'ER':
        net = ErdosRenyiGraph(N, p, directed)
        
    elif netmodel == 'RG':
        net = RandomGraph(N, L, directed)

    # 1.2) Save network into a file
    if filetype == 'text':
        fname = filename + '_re%d.dat' %re
        savetxt(join(outpath, fname), net, fmt='%d')
    elif filetype == 'binary':
        fname = filename + '_re%d.npy' %re
        save(join(outpath,fname), net)
    elif filetype == 'pajek':
        fname = filename + '_re%d.net' %re
        Save2Pajek(join(outpath,fname), net, directed=directed)


print 'Finished.'
if filetype == 'text':
    print 'Adjacency matrices saved as text files (.dat)'
elif filetype == 'binary':
    print 'Adjacency matrices saved as binary files (.npy)'
elif filetype == 'pajek':
    print 'Networks saved as Pajek-readable files (.net)'
print 'Networks saved in %s' %outpath


        
