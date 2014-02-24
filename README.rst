GAlib -- Graph Analysis library in Python / NumPy
=================================================

GAlib is a library for the analysis of graphs and complex networks in Python/NumPy. The library is very easy to install, use, modify and extend.

It's main features are...

- extensive use of NumPy to boost performance over pure Python code,
- networks are represented by their adjacency matrices (rank-2 ndarrays),
- transparent and flexible,
- includes "*Utility Scripts*" useful also for non-(Python)-programmers,
- open-source license,

... and some limitations:

- poor management of very large networks due to NumPy dependence,
- no included graph visualization tools.


INSTALLATION
------------

GAlib does not require installation but requires to have NumPy installed. The library can be downloaded from its GitHub repository:

https://github.com/gorkazl/pyGAlib

At the bottom of the right-hand menu click on the "*Download ZIP*" button. Unzip the file and copy the modules into the standard "*Site-Packages*" folder of the current Python distribution. GAlib is composed of three files: *galib.py*, *gamodels.py* and *gatools.py*. Import them as usual to start using them, e.g.: ::

>>> from galib import*
>>> from gamodels import*
>>> from gatools import*

For more details see the User Guide in the */Docs* folder.


HOW TO FIND FURTHER DOCUMENTATION
---------------------------------

While working in an interactive session and having imported the modules, use the ``help()`` function to obtain information of the library files: ::

>>> help(modulename)

This will show a list of functions included in the module. For further details regarding each function, type: ::

>>> help(functionname)

For IPython users, the help command is replaced by a question mark after the module's or function's name, e.g.: ::

>>> modulename?
>>> functionname?


CONTACT INFORMATION
-------------------

GAlib will have soon its own webpage and domain. For the moment, contact me at <galibrary@yahoo.com> for any questions, bug reports, etc.


FUTURE DEVELOPMENTS
-------------------

The current version includes core functionalities for graph analysis. Future releases will include further functionalities, for example:

* accelaration of slowest functions using the Numba package,
* further measures for weighted and/or directed networks,
* further classical (textbook) graph models,
* roles of nodes in networks with modular organization,
* extended data conversions between graph formats,
* support for sparse matrices to increase the size of networks handled, and
* ... any extensions/functions you want to include.

**Third party collaborations to extend GAlib are highly welcome!**


HISTORY OF CHANGES
------------------


February 24th, 2014
^^^^^^^^^^^^^^^^^^^

Since the official release on September 10th 2013 I have performed a few changes, some of them thanks to the feedback from a few colleagues. Thanks in particular to Miguel Lechón and to Nikos Kouvaris. Here the list of changes:

- Issues with recurrent imports solved. Only absolute imports are allowed in the modules.
- ``Degree()`` function in *galib*.py module modified to exploit properties of boolean ndarrays.
- Functions to compute roles of nodes in modular networks included to *galib.py*.
- ``BarabasiAlbert()`` function in *gamodels.py* is now always initialized with a fully connected subgraph of ``m+1`` nodes. Otherwise some hubs remained disconnected.
- ``Reciprocity()`` function in *galib.py* is now faster using boolean ndarrays. The parameter ``weighted`` has been omitted for useless and confusing.
- ``RewireNetwork()`` in *gamodels.py* has been corrected. In the very particular case of undirected graphs with assymetric link weights, weigths were not conserved. Now all nodes conserve their input intensity also in that case.
- A new module has been included: *galib_numba.py*. This is intended for the slowest functions of GAlib to be accelerated using the Numba package. Users with Numba installed can call those faster functions independently of the main galib import. For the moment I only included my main priority, a fast function for the Floyd-Warshall algorithm, ``FloydWarshall_Numba()``.





