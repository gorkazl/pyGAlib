"""
This script is meant to easily generate ensembles of scale-free networks. 
After setting of some 
parameters to personalize the results, the script generates the networks
and saves them into files for later use. The main choices are:

NETWORK PARAMETERS
- N : Size of the network (number of nodes).
- dens : The density of links of the network generated.
- gamma : The exponent of the tail, p(k) ~ k^(-gamma)
- directed : True to generate directed networks, False for undirected. In
case of directed networks, the input and output degrees are correlated, 
i.e., the input hubs are also the output hubs.

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

> cd ... /GAlib/HelperScripts

where ... is the path in which GAlib has been copied to. Then enter the
following command:

> python GenerateRandomNetworks.py

"""

from os.path import join
from numpy import save, savetxt, arange
#from gatools import Save2Pajek
from gamodels import ScaleFreeGraph


### PARAMETERS TO BE MODIFIED BY THE USER ##################################

### NETWORK OPTIONS
# Size of the network, number of nodes
N = 1000
# The density of links
dens = 0.05
# The exponent of the SF tail. Recommended between 2 and 3
gamma = 3
# Directedness. 'True' to generate directed networks, 'False' for undirected
directed = True

# Number of realizations
realiz = 10


#### OUTPUT OPTIONS
# Path where networks will be saved. The folder must already exist.
outpath = '/Work/Data/temp/'

# Base name for the files. Final name will be 'basename_reX', where X is
# the realization number
filename = 'SFNetwork'

# Type of file to be saved.
# - 'text' saves adjacency matrices as text files with .dat extension.
# - 'binary' saves adjacency matrices in numpy binary format with .npy extension.
# - 'pajek' saves networks into Pajek-readable text files with extension .net.
filetype = 'text'



############################################################################
### UNLESS CONFIDENT, DO NOT MODIFY THE CODE FROM THIS LINE ################
# 0) GIVE SOME FEEDBACK TO THE USER
print '- Number of nodes:', N
print '- Link density:', dens

if directed: print '- Directed'
else: print '- Undirected'


# 1) DO THE WORK
print '\nGenerating and saving %d scale-free networks...' %realiz

for re in xrange(int(realiz)):
    # 1.0) Feedback the state of the calculation
    if realiz <= 200:
        if re in arange(25,realiz,25):
            print re
    else:
        if re in arange(100,realiz,100):
            print re

    # 1.1) Create the network
    net = ScaleFreeGraph(N, dens, gamma, directed)

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


        
