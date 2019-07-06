# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Graph Analysis Library
======================

A package for the analysis of graphs and complex networks in Python.
Compatible with Python 2.7 and 3.X

GAlib treats networks as adjacency matrices, represented as 2D NumPy arrays,
thus taking advantage of the faster operations over pure Pythen data types.
The core of the library consists of three modules:

metrics
    Basic graph metrics (degrees, clustering, graph distance, etc)
models
    Generation of synthetic networks and randomization.
tools
    Miscelaneous helper functions.

The rest of modules provide additional functionalities:

metrics_numba
    Uses package Numba to accelerate calculation of some metrics.
models_numba
    Uses package Numba to accelerate generation of some graph models.
extra.py
    Additional measures and functionalities related to network analysis.

USING pyGAlib
^^^^^^^^^^^^^^
The library is organised in the following modules:

- *metrics.py*: Common graph metrics (degrees, clustering, graph distance, etc)
- *models.py*: Generation of synthetic networks and randomization.
- *tools.py*: Miscelaneous helper functions.
- *metrics_numba.py*: Uses the Numba package to accelerate calculation of some
metrics.
- *models_numba.py*: Uses the Numba package to accelerate generation of some
graph models.
- *extra.py*: Additional measures and functionalities related to network analysis.

Getting started
***************
Since pyGAlib depends on NumPy, it is recommended to import NumPy first.
Although this is not necessary for loading pyGAlib, NumPy functionalities and
array manipulation will be often needed. Try importing pyGAlib: ::

    >>> import numpy as np
    >>> import galib

..note::
    Importing galib imports also all functions in module *metrics.py*
    into its namespace. The rest of modules are imported separately. Therefore,
    if the import is relative those functions can be called as, e.g.,

::

    >>> import galib
    >>> ...
    >>> deg = galib.Degree(net)
    >>> C, Cnodes = galib.Clustering(net)

See that we did not have to call ``galib.metrics.Degree(net)``. In the case of
an absolute import (using an asterisk ``*``) all functions in *metrics.py* are
imported to the base namespace:  ::

    >>> from galib import *
    >>> ...
    >>> deg = Degree(net)
    >>> C, Cnodes = Clustering(net)

Example
*******
Let's generate a random graph following the Erdos-Renyi model, G(N,p), with
*N = 100* nodes and link probability *p = 0.1*:  ::

    >>> import galib
    >>> import galib.models
    >>> N = 100; p = 0.1
    >>> net = galib.models.ErdosRenyiGraph(N, p, directed=False)

Here, *net* is the adjacency matrix of the random graph represented as a
2-dimensional NumPy array. Let's calculate some basic properties.  ::

    >>> galib.Density(net)
    0.09838383838383838

As expected, the density of an Erdos-Renyi random graph is close to the
*p = 0.1* value given. We now calculate the degree of every node:  ::

    >>> deg = galib.Degree(net)
    >>> deg
    array([10,  7, 10, 10, 11,  7,  5, 11, 13, 12, 14, 13,  8, 10,  9,  8,  7,
       10, 11,  9, 11, 11,  8, 10,  5,  9, 13, 10, 13, 12, 12, 11, 11,  7,
       13, 11,  7, 10, 10,  6, 12, 10,  6, 10,  7,  6,  9, 10,  9,  9,  7,
        9,  8, 13, 10,  9,  7,  7, 11,  8, 13,  6,  7, 12, 14,  6,  5, 11,
        5, 12, 14, 14, 13,  8,  7, 12,  4, 19,  9, 13,  7, 10, 15, 15,  4,
        9,  7, 12,  7,  8, 12,  4, 11, 12,  6, 13,  6, 12, 16, 12])

The degree is returned as a numpy array of rank 1, of integer type. The function
``Clustering()`` returns both the global clustering coefficient and the local
clustering of every node:  ::

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

We compute the pair-wise graph distance between all the nodes using the
Floyd-Warshall algorithm, which is returned as a matrix ``dij`` (numpy array of
rank 2):  ::

	>>> dij = galib.FloydWarshall(net)
	>>> avlen = (dij.sum() - dij.trace()) / (N*(N-1))
	>>> avlen
	2.248080808080808

.. note::
    For calculating the distance matrix of larger networks please use
    the version of the function located in the module *metrics_numba.py*. Here,
    the function ``FloydWarshall_Numba()`` works the same but makes use of the
    Numba library to significantly speed the calculation.

Most network generators and graph metrics in pyGAlib work with directed graphs
as well. Check for the optional parameter ``directed``. Following the example
above, we generate a directed Erdos-Renyi graph and calculate its input and
output degrees for every node:  ::

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

Data I/O
********
Since GAlib is based on NumPy arrays, saving and reading of adjacency matrices,
as well as any other output of GAlib functions, can be performed using the usual
data I/O functionalities of NumPy. See for example the documentation for
functions: ``loadtxt()``, ``savetxt()``, ``load()``, ``save()`` and ``savez()``.
The *tools.py* module in pyGAlib provides also some data conversion
functionalities.

How to find further documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While working in an interactive session, after importing a module, the built-in
``help()`` function will show further details:  ::

    >>> help(modulename)

The help for galib (``help(galib)``) shows the general summary of the package
and a list of all the modules in the library. The help for each module,
e.g., ``help(galib.metrics)`` or ``help(galib.models)`` will display module
specific information and a list of all the functions in the module. For further
details regarding each function, type:  ::

    >>> help(galib.modulename.functionname)

For IPython and Jupyter notebook users the help command is replaced by a
question mark after the module's or function's name, e.g.:  ::

    >>> modulename?
    >>> functionname?

For questions, bug reports, etc, please write to <galib@Zamora-Lopez.xyz>, or
open an issue in GitHub.

License
-------

Copyright (c) 2013 - 2019, Gorka Zamora-López <gorka@Zamora-Lopez.xyz>

Released under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

or see the LICENSE.txt file.

.. note::
    Please, use the logos provided in the Branding/ folder whenever possible.

"""
from __future__ import absolute_import

from . import metrics
from .metrics import*


__author__ = "Gorka Zamora-Lopez"
__email__ = "galib@Zamora-Lopez.xyz"
__copyright__ = "Copyright 2013-2019"
__license__ = "Apache License version 2.0"
__version__ = "1.0.4"
__update__="06/07/2019"



#
