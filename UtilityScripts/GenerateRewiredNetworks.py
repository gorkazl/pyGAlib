"""
This script is meant to easily generate ensembles of randomised (rewired) 
networks that conserve the degree sequence(s) of a given network, even for
users with no or little programming experience. Users need to point the
the original network from which rewired versions are desired and they need
to specify the number of desired realizations. The program authomatically
finds whether the network is directed or undirectednot and will generate
rewired networks accordingly and save the results into files.

INPUT PARAMETERS
Input networks can be entered in three file formats: 1) a text file
containing the adjacency matrix (typically .txt or .dat), 2) a numpy binary
file containing the adjacency matrix (.npy) or 3) a text file with the 
network in Pajek format. The input parameters are:
- inputpath : Path where the input network file is located.
- inputfile : The name and extension of the file with the network.

REWIRING PARAMETERS
- realiz : Number of rewired networks to be generated.
- prew : Rewiring 'probability'. Specifies the number of iterations of the
Markov chain. Networks will be rewired prew*L times, where L is the number 
of links in the network.
- checkresult : if True, the programm will check for every rewired network
whether the degree sequence(s) have been conserved, as-well-as for the
lack of self-loops and multiple edges. If False, the programm will skip
the check and therefore run faster.

Due to the nature of the Markov chain process, prew = 2 is the minimally
acceptable value to achieve a network with all internal topology randomized,
apart from the degrees of the nodes. If computational time is not a strong
limitation, use prew = 10. Otherwise, prew = 5 is the smallest recommended
value.

OUTPUT PARAMETERS
Equivalent to the input file options, the resulting network can also be
saved into three types of files: text files, numpy binary or Pajek-readable
files. The output parameters to specify are:
- outpath : Path where the resulting networks will be saved.
- filename : Basename for the output files. Final name will be
basename_reX.dat where X is the realization number.
- filetype : Choice to save networks. 'text' saves the adjacency matrices
as text files with a .dat extension. 'binary' saves the adjacency matrices
in the numpy binary format with a .npy extension. 'pajek' saves the networks
into text files with the Pajek format and extension .net.

RUNNING THE SCRIPT
After personalizing the script and saving the changes, open a terminal
window and move to the folder containing the script, e.g.:

> cd ~/GraphAnalysisLibrary/UtilityScripts

then type the following command:

> python GenerateRewiredNetworks.py

"""

__author__ = "Gorka Zamora-Lopez" 
__email__ = "Gorka.zamora@ymail.com"
__copyright__ = "Copyright 2013-2015"
__license__ = "GPL"
__update__="07/01/2015"

from os.path import join
from numpy import*
from gatools import Save2Pajek, LoadFromPajek
from galib import Degree, Reciprocity
from gamodels import RewireNetwork

### PARAMETERS TO BE MODIFIED BY THE USER ##################################
### NETWORK OPTIONS
# Path to the original network to be rewired
inputfile = '/yourpath/GAlib/Examples/Data/Cat53_cortex.txt'

# Number of realizations
realiz = 100

# Rewiring 'probability'. If computational time is not a limitation, 
# prew = 10 is recommended. Otherwise, prew = 5 is the smallest recommended
# value.
prew = 10

# Check results. If True the program checks whether every rewired network
# conserves the degree sequence. Otherwise type False and will run faster.
checkresult = True


#### OUTPUT OPTIONS
# Path where networks will be saved. The folder must already exist.
outpath = '/yourpath/temp/'

# Base name for the files. Final name will be 'basename_PrewX_reY', where 
# X is value of prew and Y is the realization number.
filename = 'RewiredNetwork'

# Type of file to be saved.
# - 'text' saves adjacency matrices as text files with .dat extension.
# - 'binary' saves adjacency matrices in numpy binary format with .npy extension.
# - 'pajek' saves networks into Pajek-readable text files with extension .net.
filetype = 'binary'



############################################################################
### UNLESS CONFIDENT, PLEASE DO NOT MODIFY THE CODE FROM THIS LINE.
# 0) READ THE ORIGINAL NETWORK AND PREPARE BASIC PARAMETERS
name,extension = inputfile.split('.')
if extension == 'net':
    net = LoadFromPajek(inputfile)
elif extension == 'npy':
    net = load(join(inputfile)
else:
    net = loadtxt(inputfile, uint8)
N = len(net)


# 0.1) Check whether the network is directed or undirected
recip = Reciprocity(net)
if recip == 1:
    directed = False
    deg = Degree(net)
    L = 0.5 * deg.sum()
else:
    directed = True
    ink, outk = Degree(net, True)
    L = ink.sum()

# 0.1) Give some feedback to the user
if directed:
    print '- Network is directed'
else:
    print '- Network is undirected'


# 1) DO THE WORK
print '\nGenerating and saving %d rewired networks...' %realiz

for re in xrange(realiz):
    # 1.0) Feedback the state of the calculation
    if realiz <= 200:
        if re in arange(25,realiz,25):
            print re
    else:
        if re in arange(100,realiz,100):
            print re
            
    # 1.1) Rewire the network
    rewnet = RewireNetwork(net, prew, directed)

    # 1.2) Check that the rewired network is correct
    if checkresult:
        if directed:
            # Check that degrees were conserved
            newink, newoutk = Degree(rewnet, True)
            if (newoutk-outk).any(): print re, 'outk not conserved!!'
            if (ink - newink).any(): print re, 'ink not conserved!!'

            # Check that no self-loops nor multiple links were introduced
            sl = diagonal(rewnet).sum()
            Lnew = newoutk.sum()
            ml = abs(L - Lnew)
            failures = sl + ml
            if failures > 0:
                print re, 'self-loops: %d\tmultiple links: %d' %(sl,ml)

        else:
            # Check that degrees were conserved
            newdeg = Degree(rewnet)
            if (newdeg-deg).any(): print re, 'degree not conserved!!'

            # Check that no self-loops nor multiple links were introduced
            sl = diagonal(rewnet).sum()
            Lnew = 0.5 * newdeg.sum()
            ml = L - Lnew
            failures = sl + ml
            if failures > 0:
                print re, 'self-loops: %d\tmultiple links: %d' %(sl,ml)
                #break
            

    # 1.3) SAVE THE REWIRED NETWORK INTO A FILE
    if filetype == 'text':
        fname = filename + '_Prew%d_re%d.dat' %(prew, re)
        savetxt(join(outpath, fname), rewnet, fmt='%d')
    elif filetype == 'binary':
        fname = filename + '_Prew%d_re%d.npy' %(prew, re)
        save(join(outpath,fname), rewnet)
    elif filetype == 'pajek':
        fname = filename + '_Prew%d_re%d.net' %(prew, re)
        Save2Pajek(join(outpath,fname), rewnet, directed=directed)


print 'Finished.'
if filetype == 'text':
    print 'Adjacency matrices saved as text files (.dat)'
elif filetype == 'binary':
    print 'Adjacency matrices saved as binary files (.npy)'
elif filetype == 'pajek':
    print 'Networks saved as Pajek-readable files (.net)'
print 'Networks saved in %s' %outpath


