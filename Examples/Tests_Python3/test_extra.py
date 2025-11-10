# -*- coding: utf-8 -*-
# Copyright (c) 2013 - 2022, Gorka Zamora-López <galib@Zamora-Lopez.xyz>
#
# Released under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""In this file I will test all functions in module galib.extra.py.
"""

# Standard library imports
from timeit import default_timer as timer
# Third party imports
import matplotlib.pyplot as plt
from numpy import*
import numpy.linalg
import numpy.random
# Personal libraries
from galib.tools import LoadFromPajek, Save2Pajek
from galib.extra import*


##################################################################
# 0) READ THE DATA
dataroot = '../Data/'
net, labs = LoadFromPajek(dataroot + 'Cat53_cortex.net', True)
netsym = 0.5*(net+net.T)
N = len(net)


time1 = timer()
# 1) ESTIMATION OF EXPECTED CROSS-CORRELATION
# Topological similarity
fcnet = TopologicalSimilarity(net, 2.5)

plt.figure()
plt.imshow(fcnet)
plt.clim(0,1)
plt.colorbar()

# Exponential mapping
fcnet = ExponentialMapping(net, 2.5)

plt.figure()
plt.imshow(fcnet)
plt.clim(0,1)
plt.colorbar()

# Linear gaussian estimate
evs = numpy.linalg.eigvals(netsym)
evmax = evs.max()
normnet = netsym / evmax
covmat = CovarianceLinearGaussian(normnet, 0.8)
corrmat = abs(Covmat2Corrmat(covmat))

plt.figure()
plt.imshow(corrmat)
plt.clim(0,1)
plt.colorbar()


# 2) FUNCTIONAL COMPLEXITY MEASURES
# 2.1) Functional complexity
fcomp = FunctionalComplexity(corrmat)
print( 'Functional Complexity:', fcomp )

# 2.2) Tononi's functional complexity (neural complexity, sampled!!)
ncomp = NeuralComplexity_Sampled(corrmat)
print( 'Neural Complexity (sampled):', ncomp )

# 2.3) Tononi's functional complexity (neural complexity)
ncomp = NeuralComplexity(corrmat[:15,:15])
print( 'Neural Complexity:', ncomp )


time2 = timer()
print( time2 - time1, 'Seconds' )

plt.show()




#
