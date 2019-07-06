# pyGAlib – Graph Analysis Library in Python / NumPy


pyGAlib is a library to generate and study graphs and complex networks in Python. It treats networks as adjacency matrices in order to take advantage of faster NumPy
array manipulations. The library is easy to install, use, modify and extend.

Main features...

- networks represented as adjacency matrices (rank-2 ndarrays);
- extensive use of NumPy for high performance over pure Python code;
- Some functions further boosted using Numba;
- transparent and flexible: find which part of the code is doing what;
- Python 2.7 and 3.X compatible.

... and some limitations:

- management of large networks due to NumPy dependence,
- no graph visualization tools.


### INSTALLATION

Installation of pyGAlib is simple. An existing python distribution and the [pip](https://github.com/pypa/pip) package manager need to be installed. If Python was installed via the [Canopy](https://www.enthought.com/product/canopy/) or the [Anaconda](https://www.anaconda.com) distributions, then `pip` is surely installed. To check, open a terminal and type:

	$ pip --help

**- The simple installation**: pyGAlib is registered in [PyPI](https://pypi.org/project/galib/) (the Python Packages Index), therefore installation from the terminal is straight forward. `pip` will automatically take care of the  dependencies (see the *requirements.txt* file). Simply type:

	$ pip install galib 

To confirm the installation open an interactive session and try to import the library by typing `import galib`.

> **NOTE:** If you are using Python 2 and Python 3 environments, pyGAlib needs to be installed in each of the environments separately.


**- Installation from GitHub (full download)**: Visit the GitHub repository of pyGAlib [https://github.com/gorkazl/pyGAlib/](https://github.com/gorkazl/pyGAlib/) and click on the "Clone or download" button at the right hand side (the green button). Select "Download ZIP". Unzip the file, open a terminal and move to the folder, e.g.,

	$ cd ~/Downloads/pyGAlib-master/

Once on the folder that contains the *setup.py* file, type the following

	$ pip install .

Do not forget the "." at the end which means "*look for the setup.py file in the current directory*." This will check for the dependencies and install pyGAlib. To confirm installation, try running one of the  test scripts in the *Examples/* folder, e.g.,

	$ cd Examples/Test_Python2/
	$ python test_metrics_py2.py

> **NOTE**: After installation the current folder "*~/Downloads/pyGAlib-master/*" can be safely deleted, or moved somewhere else if you want to conserve the examples and the tests.

**- Installation from GitHub (lazy version)**: If [git](https://git-scm.com) is also installed in your computer, then open a terminal and type:

	$ pip install git+https://github.com/gorkazl/pyGAlib.git@master

This will install the package, that is, the content in the folder *galib/*. Other files (Examples/, README.md, LICENSE.txt, etc.) need to be downloaded manually, if wanted.



### HOW TO USE pyGAlib

The library is organised in the following modules: 

- *metrics.py*: Common graph metrics (degrees, clustering, graph distance, etc)
- *models.py*: Generation of synthetic networks and randomization.
- *tools.py*: Miscelaneous helper functions.
- *metrics_numba.py*: Uses the Numba package to accelerate calculation of some metrics.
- *models_numba.py*: Uses the Numba package to accelerate generation of some graph models.
- *extra.py*: Additional measures and functionalities related to network analysis.

#### Getting started 

Since pyGAlib depends on NumPy, it is recommended to import NumPy first. Although
this is not necessary for loading pyGAlib, NumPy functionalities and array
manipulation will be often needed. Try importing pyGAlib:

	>>> import numpy as np
	>>> import galib

> **NOTE**: Importing galib imports also all functions in module *metrics.py* into its namespace. The rest of modules are imported separately. Therefore, if the import is relative those functions can be called as, e.g., 
 
	>>> import galib
	>>> ... 
	>>> deg = galib.Degree(net)
	>>> C, Cnodes = galib.Clustering(net)
	
> See that we did not have to call `galib.metrics.Degree(net)`. In the case of an absolute import (using an asterisk `*`) all functions in *metrics.py* are imported to the base namespace:

	>>> from galib import *
	>>> ... 
	>>> deg = Degree(net)
	>>> C, Cnodes = Clustering(net)

##### Example

Let's generate a random graph following the Erdos-Renyi model, G(N,p), with
*N = 100* nodes and link probability *p = 0.1*:

	>>> import galib
	>>> import galib.models
	>>> N = 100; p = 0.1
	>>> net = galib.models.ErdosRenyiGraph(N, p, directed=False)

Here, *net* is the adjacency matrix of the random graph represented as a
2-dimensional NumPy array. Let's calculate some basic properties.

	>>> galib.Density(net)
	0.09838383838383838

As expected, the density of an Erdos-Renyi random graph is close to the *p = 0.1* value given. We now calculate the degree of every node:

	>>> deg = galib.Degree(net)
	>>> deg
	array([10,  7, 10, 10, 11,  7,  5, 11, 13, 12, 14, 13,  8, 10,  9,  8,  7,
       10, 11,  9, 11, 11,  8, 10,  5,  9, 13, 10, 13, 12, 12, 11, 11,  7,
       13, 11,  7, 10, 10,  6, 12, 10,  6, 10,  7,  6,  9, 10,  9,  9,  7,
        9,  8, 13, 10,  9,  7,  7, 11,  8, 13,  6,  7, 12, 14,  6,  5, 11,
        5, 12, 14, 14, 13,  8,  7, 12,  4, 19,  9, 13,  7, 10, 15, 15,  4,
        9,  7, 12,  7,  8, 12,  4, 11, 12,  6, 13,  6, 12, 16, 12])

The degree is returned as a numpy array of rank 1, of integer type. The function *Clustering* returns both the global clustering coefficient and the local clustering of every node:

	>>> C, Cnodes = galib.Clustering(net)
	>>> C
	0.096051227321238
	>>> Cnodes
	array([0.13333333, 0.04761905, 0.04444444, 0.08888889, 0.10909091,
       0.19047619, 0.3       , 0.12727273, 0.08974359, 0.06060606,
      	... ... ...
      	... ... ...
       0.04545455, 0.        , 0.10909091, 0.13636364, 0.        ,
       0.07692308, 0.13333333, 0.09090909, 0.08333333, 0.10606061])

We compute the pair-wise graph distance between all the nodes using the Floyd-Warshall algorithm, which is returned as a matrix *dij* (numpy array of rank 2):

	>>> dij = galib.FloydWarshall(net)
	>>> avlen = (dij.sum() - dij.trace()) / (N*(N-1))
	>>> avlen
	2.248080808080808

> **NOTE:** For calculating the distance matrix of larger networks please use the version of the function located in the module *metrics_numba.py*. Here, the function `FloydWarshall_Numba()` works the same but makes use of the Numba library to significantly speed the calculation.

Most network generators and graph metrics in pyGAlib work with directed graphs as well. Check for the optional parameter `directed`. Following the example above, we generate a directed Erdos-Renyi graph and calculate its input and output degrees for every node:

	>>> net = galib.models.ErdosRenyiGraph(N, p, directed=True)
	>>> galib.Density(net)
	0.10272727272727272
	>>> indeg, outdeg = galib.Degree(net, directed=True)
	>>> indeg
	array([17,  7,  9,  8, 11, 10,  9,  8, 13, 13,  5,  9, 13,  9, 10, 10, 13,
       10,  9,  9,  7, 11, 13, 10,  4, 15, 11, 11, 10,  6,  6,  8,  8,  8,
       11,  8,  4, 12,  8, 13, 13, 14, 12,  5,  6,  5, 16, 12,  5, 10,  9,
       13,  8,  9,  7,  8, 13, 14,  9, 18,  7, 11,  5,  4, 12,  8,  8, 10,
        7,  9, 15, 12, 14,  9, 15, 11, 13, 12, 15, 10, 11, 11, 15,  7, 10,
       13,  7, 14,  9, 16, 11, 11,  6, 18,  7,  4, 14, 12, 12, 10])
	>>> outdeg
	array([ 9, 10,  7,  9, 12,  9, 19,  9, 11, 16, 11, 12, 11, 15, 11,  6,  9,
        8, 11, 12,  9, 13,  9,  8, 11,  6,  7, 11, 11, 12, 10,  8, 11, 12,
       10, 12, 13,  8, 18, 11,  8, 13, 10,  8, 10, 10, 11,  8, 11, 11, 11,
       10, 11, 10,  9, 12,  6, 10,  7, 10, 10, 11, 15, 12, 11,  7, 10,  8,
        5, 11,  7, 11, 13,  8,  5,  6, 13, 11, 10, 13,  7,  6, 13, 11,  8,
        8, 10,  6, 10,  9, 12, 15, 11,  9, 15, 11,  7,  8, 11, 10])

##### Data I/O

Since GAlib is based on NumPy arrays, saving and reading of adjacency matrices,
as well as any other output of GAlib functions, can be performed using the
usual data I/O functionalities of NumPy. See for example the documentation for functions: `loadtxt()`, `savetxt()`, `load()`, `save()` and `savez()`. The *tools.py* module in pyGAlib provides also some data conversion functionalities.



#### Finding further documentation

While working in an interactive session, after importing a module, the built-in `help()` function will show further details:

	>>> help(modulename)

The help for galib (`help(galib)`) shows the general summary of the package and a list of all the modules in the library. The help for each module, e.g., `help(galib.metrics)` or `help(galib.models)` will display module specific information and a list of all the functions in the module.
For further details regarding each function, type:

	>>> help(galib.modulename.functionname)

For IPython and Jupyter notebook users the help command is replaced by a question mark after the module's or function's name, e.g.:

	>>> modulename?
	>>> functionname?

For questions, bug reports, etc, please write to <galib@Zamora-Lopez.xyz>, or open an issue in GitHub.


### FUTURE DEVELOPMENTS

See the TODO.md file. 
**Collaborations to extend pyGAlib are welcome.** If you have experience using *scipy.sparse*, developing community detection methods or coding graph visualization, please, please, contact me. 


### LICENSE
Copyright 2018, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

**Note:** Please, use the logos provided in the *Branding/* folder whenever possible.

-----------------------------------------------------------------
### WHAT IS NEW

##### July 6, 2019
For clarity, function `RichClub()` has been splitted into two functions: `RichClub()`and `k_Density()`. The reason is that the output of old `RichClub()` was basically the k-density for all degrees, from 0 to kmax. Now this is done by `k_Density()` and the new `RichClub()` function identifies the set of nodes (hubs) for which k-density overcomes a given value.

##### June 15, 2019
GAlib has been registered in PyPI ([https://pypi.org/project/galib/](https://pypi.org/project/galib/)). Direct installation and version management using `pip` is thus available.

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
 
See the file *CHANGES.md* for a complete history of changes.






