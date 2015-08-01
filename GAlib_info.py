"""
Graph Analysis Library
======================

A library for the analysis of graphs and complex networks in Python/NumPy.

The library contains 3 modules. See their documentation for a detailed list
of functions within each module:
- gatools.py : miscelaneous helper functions, e.g. data-type transformations.
- galib.py : basic graph descriptors, e.g. degrees and clustering.
- galib_models.py : Generation of synthetic networks and randomization.

In GAlib, graphs and networks are represented by their adjacency matrices,
that is, 2-dimensional numpy arrays. The library makes extensive use of
NumPy's array manipulation tools to compute the metrics of the networks.

For example, an empty network of N nodes is a NxN empty numpy array

>>> import numpy as np
>>> emptynet = np.zeros((N,N), dtype=uint8)

The dtype is optional and depends on whether the network will be binary or
weighted. See the numpy documentation for the available data types.

How to use GAlib
----------------

After importing GAlib you can apply its functions to an exisiting network
or generate a synthetic network for its analysis. In the following example
we create a random graph of N = 1000 nodes and perform some basic analysis.

>>> import galib
>>> import galib_models as gam
>>> N = 1000
>>> p = 0.05
>>> net = gam.ErdosRenyiGraph(N,p)

Compute the degree of every node and the mean degree using numpy a method.
>>> degree = galib.Degree(net)
>>> degree.mean()
49.603999999999999
>>>

Compute the clustering coefficient and the clustering of every node
>>> clustcoef, clustnodes = galib.Clustering(net)
>>> clustcoef
0.050222165447963818
>>> clustnodes[:10]
array([ 0.05365854,  0.05254902,  0.05325815,  0.03765227,  0.04019038,
        0.05019608,  0.0538415 ,  0.04192872,  0.04941176,  0.04542278])
>>>

GAlib works also with directed networks. Measures for directed networks
are usually specified with the optional argument 'directed'. By default,
GAlib functions come with 'directed=False'. In the following example we
generate a random directed graph and compute the input and output degrees
of its nodes.

>>> net = gam.ErdosRenyiGraph(N, p, directed=True)
>>> indegree, outdegree = galib.Degree(net, directed=True)
>>> indegree[:10]
array([48, 52, 49, 49, 55, 48, 56, 48, 48, 48])
>>> outdegree[:10]
array([61, 49, 48, 42, 51, 51, 54, 58, 52, 50])
>>>

Notice that when functions are called for directed network measures, they
often return 2-dimensional arrays that need to be adecuately unpacked as
in the example above.


Data I/O
--------

Reading networks from files and saving them into files is done using the
standard data I/O tools of numpy: loadtxt() and savetxt() are available to
work with text files. Funtions load() and save() work with the standard binary
data format of numpy. For large datasets  we recommend load() and save()
because of their better performance.

>>> savetxt('filepath.txt', net, fmt='%d)
>>> net = loadtxt('filepath.txt', dtype=uint8)

or

>>> save('filepath.npy', net)
>>> net = load('filepath.npy')

Notice that savetxt() and loadtxt() require data type formater options, while
save() and read() funcitons authomatically remember the dtype of the array.

"""