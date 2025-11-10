[![pypi version](https://img.shields.io/pypi/v/galib?logo=pypi)](https://pypi.org/project/galib/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/galib.svg?label=PyPI%20downloads)](
https://pypi.org/project/galib/)
[![Apache-2.0 License](https://img.shields.io/badge/license-Apache-blue.svg?style=flat)](http://choosealicense.com/licenses/Apache-2.0/)


# pyGAlib – Graph Analysis Library in Python / NumPy

> **Development branch**. Version 2 of pyGAlib under construction. Heavy changes expected in this branch. For download, please install the PyPI version or the one in the master branch of this repository. Both are installed via `pip`, see instructions in the corresponding README.md file. 

pyGAlib is a library to generate and study graphs and complex networks in Python. It treats networks as adjacency matrices in order to take advantage of faster NumPy
array manipulations. The library is easy to install, use, modify and extend.

Main features...

- networks represented as adjacency matrices (rank-2 ndarrays);
- extensive use of NumPy for high performance over pure Python code;
- Some functions further boosted using Numba;
- transparent and flexible: find which part of the code is doing what;
- Python 3 compatible.

... and some limitations:

- management of large networks due to NumPy dependence,
- no graph visualization tools.



### INSTALLATION

Installation of pyGAlib is simple, only the [pip](https://github.com/pypa/pip) package manager is needed. To check whether `pip` is installed, open a terminal and type:

	pip --help

> **NOTE**: If you use Anaconda (or any other third-party package manager), we recommend to install the dependencies (python>=3.6, numpy>=1.6, scipy and numba) into the target environment using Anaconda before installing pyGAlib. Otherwise, `pip` will download and install those packages directly from PyPI as well, and you won't be able to manage them through Acanconda.

#### Installing from PyPI 

pyGAlib is registered in the official *Python Package Index*, [PyPI](https://pypi.org/project/galib/) . To install, open a terminal window and type:

	python3 -m pip install galib

To confirm the installation, open an interactive session (e.g., IPython or a Notebook) and try to import the library by typing `import galib`.

#### Direct installation from GitHub 

If you have [git](https://git-scm.com) installed, you may like to install pyGAlib directly from its GitHub repository. Open a terminal and type:

	python3 -m pip install git+https://github.com/gorkazl/pyGAlib.git@master

This will only download and install the package (files in "*src/galib/*") into your current environment. If desired, additional files of the repository (e.g. the examples in the *Examples/* folder) should be downloaded manually. You can choose to install the version in another branch by replacing the '*@master*' at the end of the command by '*@branchname*' of the desired branch.

#### Installing pyGAlib in editable mode

If you want to install pyGAlib such that you can make changes to it "*on the fly*" then, visit its GitHub repository [https://github.com/gorkazl/pyGAlib/](https://github.com/gorkazl/pyGAlib/), select a branch and then click on the green "*<> Code*" button on the top right and select "Download ZIP" from the pop-up menu. Once downloaded, move the *zip* file to a target folder (e.g., "*~/Documents/myLibraries/*") and unzip the file. Open a terminal and `cd` to the resulting folder, e.g.,

	cd ~/Documents/myLibraries/pyGAlib-master/

Once on the path (make sure it contains the *pyproject.toml* file), type:

	python3 -m pip install -e .

Do not forget the "." at the end which means "*look for the pyproject.toml file in the current directory*." This will install pyGAlib such that every time changes are made to the package (located in the path chosen), these will be inmediately available. You may need to restart the IPython or Jupyter notebook session, though.



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



### HOW TO FIND FURTHER DOCUMENTATION

While working in an interactive session, after importing a module, the built-in `help()` function will show further details:

	>>> help(modulename)

The help for galib (`help(galib)`) shows the general summary of the package and a list of all the modules in the library. The help for each module, e.g., `help(galib.metrics)` or `help(galib.models)` will display module specific information and a list of all the functions in the module.
For further details regarding each function, type:

	>>> help(galib.modulename.functionname)

IPython and Jupyter notebook users, the help command is replaced by a question mark after the module's or function's name, e.g.:

	>>> modulename?
	>>> functionname?

For questions, bug reports, etc, please write to <gorka@Zamora-Lopez.xyz>, or open an issue in GitHub.


### FUTURE DEVELOPMENTS

See the *[TODO.md](https://github.com/gorkazl/pyGAlib/blob/master/TODO.md)* file in the GitHub Repository. 
**Collaborations to extend pyGAlib are welcome.** If you have experience using *scipy.sparse*, developing community detection methods or coding graph visualization, please, please, contact me. 


### LICENSE
Copyright (c) 2018, Gorka Zamora-López, <<gorka@Zamora-Lopez.xyz>>

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

##### November 10, 2025 (Release of Version 2)

Stable version 2.0 checked, validated and released.

* Python 2 support has been dropped. Only Python 3 compatibility will be developed and maintained from now on.
* The library has been reshaped to be compliant with the modern [PyPA specifications](https://packaging.python.org/en/latest/specifications/).
* [Hatch](https://hatch.pypa.io/latest/) was chosen as the tool to build and publish the package. See the *pyproject.toml* file. 
* Bug fixes to adapt to the various changes in Python and NumPy since last release.
* Sample and validation scripts in the "*Examples/*" folder revised and adapted to recent changes in Python and NumPy. 

##### March 14, 2024
Small bugs fixed:

- Normalization of `galib.metrics.Modularity()` function corrected.
- Fixed the new  aliases for `int` and `float` in *Numpy*. All arrays are now declared as `np.int64` or `np.float64`, and individual numbers as standard Python `int` or `float`. 

##### February 7, 2022
Minor bug fixes. A remaining Python 2 to Python 3 conversion error was fixed, since standard library function `range()` no longer returns a list, but an iterator object.

##### June 15, 2020
Docstrings corrected. Function `k_DensityW()` was added in module *metrics.py* to calculate the k-density in networks with weighted links, which is needed to evaluate potential formation of rich-club structures in weigthed networks.

##### July 12, 2019
Version 1.1.0 released. Section for classic and deterministic graphs added to the *models.py* module. New generators `PathGraph()`, `StarGraph()` and `CompleteGraph()` included.

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
 
See the file *CHANGELOG.md* for a complete history of changes.






