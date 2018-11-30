GAlib – Graph Analysis library in Python / NumPy
================================================

-----

**DEPRECATED BRANCH.** This is the older version of GAlib and is left here for legacy reasons only. The branch is "frozen" and no further pull requests will be accepted.
See the '*master*' branch for the library reshaped into a proper Python package.

-----

GAlib is a library for the analysis of graphs and complex networks in Python. It treats networks as adjacency matrices in order to take advantage of faster NumPy
array manipulations. The library is very easy to install, use, modify and extend.

Main features...

- extensive use of NumPy to boost performance over pure Python code;
- networks represented as adjacency matrices (rank-2 ndarrays);
- transparent and flexible: easily find which part of the code is doing what;
- "*Utility Scripts*" useful also for non-(Python)-programmers;
- open-source license.

... and some limitations:

- limited management of large networks due to NumPy dependence,
- no graph visualization tools.


### INSTALLATION

Installation of GAlib is very simple, it only requires to manually copying the modules of the library into a special folder named "site-packages" and creating a pointer for python to know where to find the library.

Requirements: Python 2.7, NumPy, SciPy, Numba (optional):

1) Locate the "site-packages" folder for the Python 2.7 distribution in which GAlib will be installed. For example, **Canopy** users will find it in:  "*/~/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/*"
**Anaconda** users with a Python 2.7 environment installed (e.g. named as 'py2') will find the folder in, for example:
"*/anaconda3/envs/py2/lib/python2.7/site-packages/.*"
2) Visit the GitHub repository for [pyGAlib](https://github.com/gorkazl/pyGAlib) and download the library clicking on the "Clone or Download" button (in green). Select "Download ZIP".
3) Unzip and remane the folder as "pyGAlib" (drop the "-master").
4) Move this folder into the "site-packages" folder previously identified.
5) Using a text editor create an empty file named *pyGAlib.pth* and save it in the "site-packages" folder. Add one line to this file containing the path to the pyGAlib library.
6)	Check installation opening an interactive python session and typing `import galib`.


### HOW TO FIND FURTHER DOCUMENTATION

While working in an interactive session and having imported the modules, use the `help()` function to obtain information of the library files:

	>>> help(modulename)

This will show a list of functions included in the module. For further details regarding each function, type:

	>>> help(functionname)

For IPython users, the help command is replaced by a question mark after the module's or function's name, e.g.:

	>>> modulename?
	>>> functionname?

Please write me to <galib@Zamora-Lopez.xyz> for any questions, bug reports, etc.


### FUTURE DEVELOPMENTS

The current version includes core functionalities for graph analysis. Future releases will include:

* accelaration of slowest functions using the Numba package,
* further measures for weighted and/or directed networks,
* further classical (textbook) graph models,
* extended data conversions between graph formats,
* support for sparse matrices to increase the size of networks handled, and
* ... any extensions/functions you want to include.

**Third party collaborations to extend GAlib are highly welcome!**


----------------------------------------------------
### HISTORY OF CHANGES

##### November 29 2018
The old version GAlib, consisting of a simple library of various modules, has been deprecated and moved into this separate branch '*old_GAlib*' for reference reasons. GAlib has been reshaped into a proper Python Package, continued in the '*master*' branch.

##### March 30th 2018
New step towards accelerating slowest functions:
- Function `MatchingIndex_Numba()` added to module *galib_numba.py*
- New module *gamodels_numba.py* has been added to GAlib. This module will contain faster network generation and rewiring functions, accelerated using Numba. The module is started with a function `RandomGraph_Numba()`, with same functionality as its equivalent in the original *gamodels.py* module but significantly faster.

##### March 29th 2018
Replaced README.rst by a README.md file and minor issues corrected in *galib.py* module.


##### May 3rd 2017
GAlib has a new module named *galib_extra.py*. This module does not provide core functionalities for the analysis of complex networks but it is intended to host measures I have developed as the result of my research, and third-party measures and functionalities which I have adapted myself into Python. These measures are all related to complex network analysis and network dynamics but they may be difficult to understand as typical graph analysis tools. It includes functions to compute the **Functional Complexity** and the **Topological Similarity** measures we have recently published. See:

- Zamora-Lopez et al. “[Functional complexity emerging from anatomical constraints in the brain: the significance of network modularity and rich-clubs](https://doi.org/10.1038/srep38424)” *Sci. Reps*. **6**:38424 (2016)
-  Bettinardi et al “[How structure sculpts function: Unveiling the contribution of anatomical connectivity to the brain’s spontaneous correlation structure](https://doi.org/10.1063/1.4980099)” *Chaos* **27**, 047409 (2017).


##### August 6th 2016
Function ``gamodels.RandomGraph()`` improved. Method to choose source and target nodes to connect changed. Originally I would use the official ``random.choice()`` function to select nodes from the list of nodes. Now it uses ``N * numpy.random.rand()`` to choose the index. The latter was slower years ago but with newer Numpy versions, it has become a faster option.

##### September 14th 2015
Function ``gatools.AllBipartitions()`` improved and extended. (i) Internal dependence with function ``gatools.AllCombinations()`` removed, (ii) function is now faster, and (iii) optional parameter *’comblength’* has been introduced in case that only the bipartitions of given length are desired.

##### August 24th 2015
Solved a bug in function ``galib.RichClub()`` that affected only the calculations with the option *’average’*.


##### August 1st 2015
- Pulled and merged into master branch changes by Schmigu (see update on June 14th).
- The function ``ModularInhomogeneousGraph()`` was added to *gamodels.py* module to generate random modular networks in which the size and/or the density of each module can be specified.
- Example in ``gatools.ExtractSubmatrix()`` corrected. The example used a 3x3 array when a 4x4 array was required.


##### June 14th 2015
- Started adding unit-tests.
- Removed ``gatools.ArrayCompare()`` since the function was broken. The functionality is covered by ``numpy.array_equal()``
- Fixed bug in ``gatools.Quartiles()``

##### April 13th 2015
Function to compute the Fisher-corrected mean correlation values, ``MeanCorrelation()``, has been added to *gatools.py* module.

##### December 23rd 2014
Following publication of Klimm et al. “[Individual node’s contribution to the mesoscale of complex networks](http://iopscience.iop.org/article/10.1088/1367-2630/16/12/125006/meta)” *New J. Phys.* **16**:125006 (2014), functions to compute the roles of nodes have been included into the main library *galib.py* Individual functions for each of the four new measures  are available:

- ``GlobalHubness()``, ``LocalHubness()``, ``NodeParticipation()``, ``NodeDispersion()``.
- Function name ``ParticipationIndex()`` changed to ``NodePArticipation()``.
- Additionally, the function ``NodeRoles()`` returns all the four measures at once.
- Functions to compute the *participation matrix* and the *participation vectors* of every node have been included.


##### February 24th 2014
Since the official release on September 10th 2013, I have performed several changes, some of them thanks to the feedback from a few colleagues. Thanks in particular to Miguel Lechón and to Nikos Kouvaris. Here the list of changes:

- Issues with recurrent imports solved. Only absolute imports are allowed in the modules.
- ``Degree()`` function in *galib*.py module modified to exploit properties of boolean ndarrays.
- Functions to compute roles of nodes in modular networks included to *galib.py*.
- ``BarabasiAlbert()`` function in *gamodels.py* is now always initialized with a fully connected subgraph of ``m+1`` nodes. Otherwise some hubs remained disconnected.
- ``Reciprocity()`` function in *galib.py* is now faster using boolean ndarrays. The parameter ``weighted`` has been omitted for useless and confusing.
- ``RewireNetwork()`` in *gamodels.py* has been corrected. In the very particular case of undirected graphs with assymetric link weights, weigths were not conserved. Now all nodes conserve their input intensity also in that case.
- A new module has been included: *galib_numba.py*. This is intended for the slowest functions of GAlib to be accelerated using the Numba package. Users with Numba installed can call those faster functions independently of the main galib import. For the moment I only included my main priority, a fast function for the Floyd-Warshall algorithm, ``FloydWarshall_Numba()``.
