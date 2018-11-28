"""
Graph Analysis Library
======================

A package for the analysis of graphs and complex networks in Python.
This is the Python 2.7 compatible version.

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
    Uses package Numba to accelerate slowest metric calculations.
models_numba
    Uses package Numba to accelerate slowest model generators.
extra.py
    Additional measures and functionalities related to network analysis.

Using GAlib
-----------

Since GAlib depends on NumPy, it is recommended to import NumPy first. Although
this is not necessary for loading GAlib, NumPy functionalities and array
manipulation will be often needed. A simple import of GAlib  ::

>>> import numpy as np
>>> import galib

will load all functions within the module metrics.py into the namespace of
galib. For example, the function to calculate the degree of every node is to
be called as ``galib.Degree(adjmatrix)``. An absolute import  ::

>>> import numpy as np
>>> import galib

loads all functions directly in the namespace of the session, thus the
previous function is simply called as ``Degree(adjmatrix)``. The rest of modules
need to be imported separately.

An empty network of N nodes is a NxN empty numpy array  ::

>>> emptynet = np.zeros((N,N), dtype=uint8)

The following example generates a Erdos-Renyia random graph, G(N,p), of
N nodes and link probability p and calculates some basic metrics:  ::

>>> import galib
>>> import galib.models as gam
>>> net = gam.ErdosRenyiGraph(N, p, directed=False)
>>> deg = galib.Degree(net)
>>> Cnet, Cnodes = galib.Clustering(net)
>>> dij = galib.FloydWarshall(net)
>>> avlen = (dij.sum() - dij.trace()) / (N*(N-1))

Here, ``net`` is the adjacency matrix of a random graph represented as a
2-dimensional NumPy array. ``deg`` is an array containing the degree of every
node, ``Cnet`` is the clustering coefficient of the network and ``Cnodes`` the
local clustering of each node. Then, ``dij`` is the pair-wise distance matrix
given as a 2-dimensional NumPy array and ``avlen`` is the average pathlength
of ``net``.

Most network generators and graph metrics work with directed graphs. Check for
the optional parameter ``directed``. In the example above, the input and output
degrees of each node would be calculated as:  ::

>>> net = gam.ErdosRenyiGraph(N, p, directed=True)
>>> indeg, outdeg = galib.Degree(net, directed=True)

Data I/O
--------

Since GAlib is based on NumPY arrays, saving and reading of adjacency matrices,
as well as any other output of GAlib functions, can be performed using the
usual data I/O functionalities of NumPy. See for example the NumPY
documentation for functions: ``loadtxt()``, ``savetxt()``, ``load()``,
``save()`` and ``savez()``.

Further information
-------------------

To see the list of all functions available within each module, use the
standard help in an interactive session, e.g.,  ::

>>> import galib
>>> help(galib.metrics)
>>> import galib.models as gam
>>> help(gam)

Same, to find further details of every function within each module:, e.g.,  ::

>>> help(galib.metrics.Degree)
>>> help(gam.ScaleFreeGraph)

In the IPython interactive session or in a Jupyter Notebook, help is requested
by typing the module of function name, and adding an interrogation mark, e.g., ::

>>> galib.metrics?
>>> gam.ScaleFreeGraph?

License
-------

See LICENSE.txt file.
Copyright (C) 2013 - 2018, Gorka Zamora Lopez, Ph.D.

"""
from __future__ import absolute_import

from . import metrics
from .metrics import*


# __all__ = ['metrics']
# __all__ = ['metrics.Degree', 'metrics.Clustering']


#
