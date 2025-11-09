### HISTORY OF CHANGES

##### March 14, 2024 (version 1.1.4)
Small bugs fixed:

- Normalization of `galib.metrics.Modularity()` function corrected.
- Fixed the new  aliases for `int` and `float` in *Numpy*. All arrays are now declared as `np.int64` or `np.float64`, and individual numbers as standard Python `int` or `float`. 

##### February 7, 2022
Minor bug fixes. A remaining Python 2 to Python 3 conversion error was fixed, since standard library function `range()` no longer returns a list, but an iterator object.

##### June 15, 2020
Docstrings corrected. Function `k_DensityW()` was added in module *metrics.py* to calculate the k-density in networks with weighted links, which is needed to evaluate potential formation of rich-club structures in weigthed networks.

##### July 12, 2019
Version 1.1.0 released. Section for classic and deterministic graphs added to the *models.py* module. New generators `PathGraph()`, `StarGraph()` and `CompletGraph()` included.

##### July 6, 2019
For clarity, function `RichClub()` has been splitted into two functions: `RichClub()`and `k_Density()`. The reason is that the output of old `RichClub()` was basically the k-density for all degrees, from 0 to kmax. Now this is done by `k_Density()` and the new `RichClub()` function identifies the set of nodes (hubs) for which k-density overcomes a given value.

##### June 15, 2019
GAlib has been registered in PYPI ([https://pypi.org/project/galib/](https://pypi.org/project/galib/)). Direct installation and version management using `pip` is thus available.

##### January 29, 2019
New in Version 1.0.1:

- Minor corrections overall.
- Function to generate Ravasz-Barabási hierarchical networks added to *models.py* module.

##### December 3, 2018
Release of Version 1.0.0 of pyGAlib. The library is now shaped as a proper Python package and is installable using standard tools. The structure of the package has been renewed and contains the following modules: 

- *metrics.py*: Common graph metrics (degrees, clustering, graph distance, etc)
- *models.py*: Generation of synthetic networks and randomization.
- *tools.py*: Miscelaneous helper functions.
- *metrics_numba.py*: Uses the Numba package to accelerate calculation of some metrics.
- *models_numba.py*: Uses the Numba package to accelerate generation of some graph models.
- *extra.py*: Additional measures and functionalities related to network analysis.


##### March 30, 2018
New step towards accelerating slowest functions:

- Function `MatchingIndex_Numba()` added to module *galib_numba.py*.
- New module *gamodels_numba.py* has been added to GAlib. This module will contain faster network generation and rewiring functions, accelerated using Numba. The module is started with a function `RandomGraph_Numba()`, with same functionality as its equivalent in the original *gamodels.py* module but significantly faster.


##### March 29, 2018
Replaced README.rst by a README.md file and minor issues corrected in *galib.py* module.


##### May 3, 2017
GAlib has a new module named *galib_extra.py*. This module does not provide core functionalities for the analysis of complex networks but it is intended to host measures I have developed as the result of my research, and third-party measures and functionalities which I have adapted myself into Python. These measures are all related to complex network analysis and network dynamics but they may be difficult to understand as typical graph analysis tools. It includes functions to compute the **Functional Complexity** and the **Topological Similarity** measures we have recently published. See:

- Zamora-Lopez et al. “[Functional complexity emerging from anatomical constraints in the brain: the significance of network modularity and rich-clubs](https://doi.org/10.1038/srep38424)” *Sci. Reps*. **6**:38424 (2016)
-  Bettinardi et al “[How structure sculpts function: Unveiling the contribution of anatomical connectivity to the brain’s spontaneous correlation structure](https://doi.org/10.1063/1.4980099)” *Chaos* **27**, 047409 (2017).


##### August 6, 2016
Function ``gamodels.RandomGraph()`` improved. Method to choose source and target nodes to connect changed. Originally I would use the official ``random.choice()`` function to select nodes from the list of nodes. Now it uses ``N * numpy.random.rand()`` to choose the index. The latter was slower years ago but with newer Numpy versions, it has become a faster option.


##### September 14, 2015
Function ``gatools.AllBipartitions()`` improved and extended. (i) Internal dependence with function ``gatools.AllCombinations()`` removed, (ii) function is now faster, and (iii) optional parameter *’comblength’* has been introduced in case that only the bipartitions of given length are desired.


##### August 24, 2015
Solved a bug in function ``galib.RichClub()`` that affected only the calculations with the option *’average’*.


##### August 1, 2015
- Pulled and merged into master branch changes by Schmigu (see update on June 14th).
- The function ``ModularInhomogeneousGraph()`` was added to *gamodels.py* module to generate random modular networks in which the size and/or the density of each module can be specified.
- Example in ``gatools.ExtractSubmatrix()`` corrected. The example used a 3x3 array when a 4x4 array was required.


##### June 14, 2015
- Started adding unit-tests.
- Removed ``gatools.ArrayCompare()`` since the function was broken. The functionality is covered by ``numpy.array_equal()``
- Fixed bug in ``gatools.Quartiles()``


##### April 13, 2015
Function to compute the Fisher-corrected mean correlation values, ``MeanCorrelation()``, has been added to *gatools.py* module.


##### December 23, 2014
Following publication of Klimm et al. “[Individual node’s contribution to the mesoscale of complex networks](http://iopscience.iop.org/article/10.1088/1367-2630/16/12/125006/meta)” *New J. Phys.* **16**:125006 (2014), functions to compute the roles of nodes have been included into the main library *galib.py* Individual functions for each of the four new measures  are available:

- ``GlobalHubness()``, ``LocalHubness()``, ``NodeParticipation()``, ``NodeDispersion()``.
- Function name ``ParticipationIndex()`` changed to ``NodePArticipation()``.
- Additionally, the function ``NodeRoles()`` returns all the four measures at once.
- Functions to compute the *participation matrix* and the *participation vectors* of every node have been included.


##### February 24, 2014
Since the official release on September 10th 2013, I have performed several changes, some of them thanks to the feedback from a few colleagues. Thanks in particular to Miguel Lechón and to Nikos Kouvaris. Here the list of changes:

- Issues with recurrent imports solved. Only absolute imports are allowed in the modules.
- ``Degree()`` function in *galib*.py module modified to exploit properties of boolean ndarrays.
- Functions to compute roles of nodes in modular networks included to *galib.py*.
- ``BarabasiAlbert()`` function in *gamodels.py* is now always initialized with a fully connected subgraph of ``m+1`` nodes. Otherwise some hubs remained disconnected.
- ``Reciprocity()`` function in *galib.py* is now faster using boolean ndarrays. The parameter ``weighted`` has been omitted for useless and confusing.
- ``RewireNetwork()`` in *gamodels.py* has been corrected. In the very particular case of undirected graphs with assymetric link weights, weigths were not conserved. Now all nodes conserve their input intensity also in that case.
- A new module has been included: *galib_numba.py*. This is intended for the slowest functions of GAlib to be accelerated using the Numba package. Users with Numba installed can call those faster functions independently of the main galib import. For the moment I only included my main priority, a fast function for the Floyd-Warshall algorithm, ``FloydWarshall_Numba()``.
