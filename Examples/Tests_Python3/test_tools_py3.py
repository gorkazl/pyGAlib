# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2019, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""In this file I will test all functions in the Python3 version of GAlib
(the module gatools.py) and make sure they work.
"""

# Standard library imports
import os, os.path
from timeit import default_timer as timer
# Third party imports
from numpy import*
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy.random
# Personal libraries
from galib import*
from galib.tools import*

##################################################################
# 1) I/O  AND DATA CONVERSIONS
currdir = os.getcwd()
currdir = os.path.split(currdir)[0]
dataroot = os.path.join(currdir, 'Data/')

time1 = timer()
net, labs = LoadFromPajek(dataroot + 'Cat53_cortex.net', getlabels=True)

# netsym = 0.5*(net + net.T)
netsym = SymmetriseMatrix(net)
print(Reciprocity(net), Reciprocity(netsym))
nzidx = nonzero(net)
net = net.astype(float32)
net[nzidx] += 0.5

Save2Pajek(dataroot + 'spam.net', net, labels=labs, directed=True, weighted=True)

labs = LoadLabels(dataroot + 'Authorlist.txt')

SaveLabels(dataroot + 'spam_labels.txt', labs)

# Define the partition
visual = arange(16)
audit = arange(16,23)
somatomotor = arange(23,39)
frontolimbic = arange(39, 53)
partition = [visual,audit,somatomotor,frontolimbic]
ncoms = len(partition)

SavePartition(dataroot + 'spam_CatPartition.txt', partition)
newpart = LoadPartition(dataroot + 'spam_CatPartition.txt')

# Slicing (adjacency) matrices
subnet1 = ExtractSubmatrix(net, partition[0])
subnet2 = ExtractSubmatrix(net, partition[0], partition[1])

# Compute the Laplacian matrix
lapnet = LaplacianMatrix(net)
print(lapnet.sum(axis=1))
print(lapnet.diagonal())

# Clean paths
Dij, bcnodes, allpaths, allcycles = PathsAllinOne(net)
newpaths = allpaths.copy()
for i in range(1,5):
    newpaths[i] += newpaths[i][:1000]
    CleanPaths(newpaths[i])
    print(len(allpaths[i]), len(newpaths[i]))


# 2) ARRAY AND MATRIX COMPARISONS
print('Hamming distance...')
spam1 = zeros(100, float)
spam2 = ones(100, float)
print( HammingDistance(spam1, spam2) )

spam2[:50] = 0
print( HammingDistance(spam1, spam2) )

spam1[10:50] = 3.5
print( HammingDistance(spam1, spam2) )


# 3) ADDITIONAL MATH
# NonZeroMin
print( 'Min values:', NonZeroMin(spam1), NonZeroMin(spam2) )

# Cumulative distribution
data = numpy.random.rand(10000)
xdata, ydata = CumulativeDistribution(data, 10, range=[0,1], normed=True, centerbins=True)

plt.figure()
plt.plot(xdata,ydata, '.')
plt.grid()

# Factorial function
for i in range(10):
    print( i, Factorial(i) )

# Binomial
for i in range(11):
    print( i, BinomialCoefficient(10,i) )

# Quartiles
q1, q2, q3 = Quartiles(data)
print( 'quartiles: %2.5f %2.5f %2.5f' %(q1,q2,q3) )

# Permuations, combinations, etc.
newdata = [1,2,4,'spam']
allperms = AllPermutations(newdata)
print( 'Permutations\n', allperms )

print('All combinations' )
for i in range(0,5):
    print( i, AllCombinations(newdata,i) )

allbips = AllBipartitions(newdata)
print( 'All bipartitions\n', allbips )

# Mean correlation
print( 'Mean correlation:', MeanCorrelation(data) )

time2 = timer()
print( time2 - time1, 'Seconds' )

plt.show()




#
