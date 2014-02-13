"""
============================
SYNTHETIC NETWORK GENERATORS
============================

This module contains functions to generate typical synthetic networks,
including random networks and methods to rewire networks.

RANDOM NETWORK GENERATORS
=========================
Lattice1D
    Returns a ring lattice where each node connects its 2z closest neighbours.
WattsStrogatzGraph
    Returns a small-world network as in the Watts & Strogatz model.
ErdosRenyiGraph
    Returns a random graphs following the Erdos & Renyi model.
RandomGraph
    Returns a random graph with N nodes and L links.
BarabasiAlbertGraph
    Returns a scale-free network after the Barabasi & Albert model.
ScaleFreeGraph
    Returns a scale-free graph of given size and exponent.

NETWORK REWIRING/RANDOMIZATION ALGORITHMS
=========================================
RewireNetwork
    Returns a network with links rewired with same degrees as net.

HIERARCHICAL AND MODULAR (HM) NETWORK MODELS
============================================
HMpartition
    Returns a partition of nodes for a hierarchical/modular random network..
HMRandomNetwork
    Returns a random hierarchical/modular network of given shape.
RavaszBarabasiModel
    Returns a hierarchical/modular network of the Ravasz & Barabasi model.
"""

import numpy as np
import numpy.random
import random

######################################################################
"""RANDOM NETWORK GENERATORS"""
def Lattice1D(N,z):
    """Returns a ring lattice where each node connects its 2z closest 
    neighbours (left and right from the node).
    
    Parameters
    ----------
    N : integer
        Size of the network (number of nodes).
    z : integer
        Diameter of the neighbourhood that every node connects with.
        For a given z, every node has degree 2*z in the resulting lattice.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the 1-dimensional lattice.
    """
    # 0) SECURITY CHECK
    assert z <= int(N/2), 'Largest possible z in N/2: %d' %(N/2)
    
    if z == 0:
        return np.zeros((N,N), np.uint8)

    # 1) CREATE THE LATTICE
    adjmatrix = np.zeros((N,N), np.int)
    
    # 1.1) Create the first row according to the number of neighbours
    adjmatrix[0,1:z+1] = 1
    adjmatrix[0,-z:] = 1
    
    # 1.2) Use numpy.roll() to copy rotated version of the first row
    for i in xrange(1,N):
        adjmatrix[i] = np.roll(adjmatrix[0],i)

    return adjmatrix

def WattsStrogatzGraph(N, z, prew, lattice=None):
    """Returns a small-world network as in the Watts & Strogatz model.
    
    See Watts & Strogatz, Nature 393(4) (1998) for details of the model.
    
    The algorithm rewires, clockwise, the links of the first neighbours, 
    of the second neighbours, etc. with probability prew. If a link is 
    chosen for rewiring, it can't create a self-loop nor a multiple edge.
    
    Parameters
    ----------
    N : integer
        Size of the network (number of nodes).
    z : integer
        Diameter of the neighbourhood that every node connects with.
        For a given z, every node has degree k = 2*z in the resulting lattice.
    prew : float, between 0 and 1.
        Probability that links of the 1D lattice to be rewired.
    lattice : ndarray of rank-2, optional.
        The adjacency matrix of a 1D lattice that has to be rewired.
        When several realizations of the model are desired, the initial
        lattice, usually the output of Lattice1D() function, can be passed
        what avoids having to create a new lattice for every realization and 
        saving computational time. Recommended
        
    Returns
    -------
    adjmatrix : ndarray of rank-2 and integer type.
        The adjacency matrix of the rewired 1-dimensional lattice.
        
    Notes
    -----
    1) Unlike the Watts & Strogatz model in which k is the total number of 
    neighbours of a node, here z represents the largest z-neighbourhood
    with which a node is connected (left and right). Hence, comparing the
    two parameters, k = 2*z.
    
    2) When several realizations of the model are desired, the initial
    lattice, usually the output of Lattice1D() function, can be passed
    to the optional argument 'lattice'. This avoids having to create a 
    new lattice for every realization and saves computational time.
    Recommended for several realizations of large networks.
    
    3) The present algorithm is only valid for graphs. In order to obtain
    directed small-world networks, use the RewireNetwork() function on a
    previously created ring-lattice.
    
    See Also
    --------
    Lattice1D : Generates a 1D lattice ring network.
    RewiredNetwork : Randomized the links of a network with equal probability.
    """
    
    # 0) SECURITY CHECKS
    assert prew >= 0 and prew <= 1, \
        "Probability 'prew' out of bounds. Insert value between 0 and 1"
    
    # 1) CREATE OR CHECK THE VALIDITY OF THE INITIAL RING-LATTICE
    if lattice==None:
        assert z <= int(N/2), 'Largest possible z in N/2: %d' %(N/2)
        adjmatrix = Lattice1D(N,z)
    else:
        Nlatt = len(lattice)
        assert N == Nlatt, 'Size N not aligned with size of given lattice.'   
        adjmatrix = lattice.copy()
    
    # 2) REWIRE THE LINKS CLOCKWISE AND BY RANK OF NEIGHBOURHOOD
    # List of nodes from which to randomly choose
    nodelist = range(N)
    # For each rank of neighbourhood...
    for k in xrange(1,z+1):
        # Choose every node in clockwise direction
        for i in xrange(N):
            # 2.1) Rewire the link with neighbour i+k with probability prew
            if numpy.random.rand() <= prew:
                done = False
                while not done:
                    target = random.choice(nodelist)
                    # Avoid self-loops and double links
                    if target == i: continue
                    if adjmatrix[i,target]: continue
                    
                    # Else, remove the old link and place the new one.
                    if i+k < N:
                        adjmatrix[i,i+k] = 0
                        adjmatrix[i+k,i] = 0
                    else:
                        adjmatrix[i,i+k-N] = 0
                        adjmatrix[i+k-N,i] = 0

                    adjmatrix[i,target] = 1
                    adjmatrix[target,i] = 1
                    done = True

    return adjmatrix

def ErdosRenyiGraph(N, p, directed=False, selfloops=False, outdtype=np.uint8):
    """Returns a random graphs following the Erdos & Renyi model.
    
    In an Erdos-Renyi graph a link happens with probability p. Therefore,
    different relizations of the graph may have different number of links.
    In the ensemble average, the expected number of links of ER graphs
    is <L> = 1/2 * p * N**2 if self-loops are permitted, or
    <L> = 1/2 * p * N * (N-1) if no self-loops are permitted.
    
    Parameters
    ----------
    N : integer
        The size of the network (number of nodes).
    p : float between 0 and 1
        Probability of links.
    directed : Boolean
        True if a directed graph is desired, False if an undirected graph is 
        desired.
    selfloops: Boolean
        True if self-loops are allowed, False otherwise.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN
        The adjacency matrix of the generated random graph.
    
    See Also
    --------
    RandomGraph : Random graphs with specified number of links
    """
    # 0) SECURITY CHECKS. Make sure that adequate parameters are given.
    assert p >= 0 and p <= 1, \
        "Probability 'p' out of bounds. Insert value between 0 and 1"

    # 1) FOR DIRECTED GRAPHS
    # 1.1) Create a random matrix with normally distributed values
    adjmatrix = numpy.random.rand(N,N)
    # 1.2 Select the entries with value <= p
    adjmatrix = np.where(adjmatrix <= p, 1, 0).astype(outdtype)

    # 2) FOR UNDIRECTED GRAPHS
    if not directed:
        adjmatrix[np.tril_indices(N)] = 0
        adjmatrix += adjmatrix.T

    # 3) Remove the diagonal if no self-loops are desired
    if not selfloops:
        adjmatrix[np.diag_indices(N)] = 0

    return adjmatrix

def RandomGraph(N, L, directed=False, selfloops=False):
    """Returns a random graph with N nodes and L links.
    
    Similar to an Erdos-Renyi (ER) graph with probability p = rho, where
    rho is the density of links. In ER graphs the total number of links
    varies in different realizations what is unsuitable to compare empirical
    networks with their random counterparts. RandomGraph() allows to create 
    random graphs that always have the same number of links. The precise 
    formula for rho depends on the options given. For an undirected graph 
    with no self-loops allowed, rho = 1/2 * L / (N*(N-1)). 

    Parameters
    ----------
    N : integer
        The size of the network (number of nodes).
    L : integer
        Number of links of the resulting random network.
    directed : Boolean
        True if a directed graph is desired. False, for an undirected graph.
    selfloops: Boolean
        True if self-loops are allowed, False otherwise.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2, size NxN and dtype = int.
        The adjacency matrix of the generated random graph.
    
    Notes
    -----
    Make sure to specify the right number of links for the type of graph
    desired. Keep in mind that directed graphs have twice as many links 
    (arcs) as graphs (edges).
        
    See Also
    --------
    ErdosRenyiGraph : Random graphs with given link probability.
    """
    # 0) SECURITY CHECKS. Make sure L is not too large.
    if directed:
        if selfloops:
            maxL = N**2
            assert L <= maxL, \
                'L out of bounds. For the options given, max(L) = N**2 = %d' %maxL
        else:
            maxL = N*(N-1)
            assert L <= maxL, \
                'L out of bounds. For the options given, max(L) = N*(N-1) = %d' %maxL
    else:
        if selfloops:
            maxL = 0.5*N*(N+1)
            assert L <= maxL, \
                'L out of bounds. For the options given, max(L) = 1/2*N*(N+1) = %d' %maxL
        else:
            maxL = 0.5*N*(N-1)
            assert L <= maxL, \
                'L out of bounds. For the options given, max(L) = 1/2*N*(N-1) = %d' %maxL

    # 1) INITIATE THE MATRIX AND HELPERS    
    adjmatrix = np.zeros((N,N), int)
    nodelist = np.arange(N)
    counter = 0

    # 2) GENERATE THE MATRIX
    while counter < L:
        # 2.1) Pick up two nodes at random
        source = random.choice(nodelist)
        target = random.choice(nodelist)

        # 2.2) Check if they can be linked, otherwise look for another pair
        if adjmatrix[source,target] == 1: continue
        if not selfloops and source == target: continue

        # 2.3) If the nodes are linkable, place the link
        adjmatrix[source,target] = 1
        if not directed:
            adjmatrix[target,source] = 1
        
        counter += 1

    return adjmatrix

def BarabasiAlbertGraph(N, m0, m, outdtype=np.uint8):
    """Returns a scale-free network after the Barabasi & Albert model.
    
    The Barabasi and Albert model (Science 286, 1999) creates networks with
    a scale-free degree distribution of exponent gamma = -3 by a growth 
    process with preferential attachment. Given an initial small connected 
    network of m0 nodes, at each iteration a new node is included that 
    connects to the existing nodes with probability proportional to their
    degree. 
    The configuration of the initial network is arbitrary. Here, we use
    a random graph with mean-degree <k> = m as a seed graph.
    
    Parameters
    ----------
    N : integer
        Number of nodes of the final network.
    m0 : integer
        Size of the seed random graph that initiates the network.
        Must be 2 or larger.
    m : integer
        Number of links that each new node makes during the growing process.
        The condition m <= m0 must hold.
    outdtype : numpy compatible data type, optional.
        The data type of the resulting adjacency matrix of the scale-free
        network.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN
        The adjacency matrix of the generated scale-free network.
    """
    
    # 0) SECURITY CHECKS
    assert m0 >= 2, 'Value not accepted. m0 must be 2 or larger.'
    assert m <= m0, 'Value not accepted. m can only take values m <= m0.'

    # 1) INITIATE THE NETWORK AS A RANDOM GRAPH OF m0 NODES AND MEAN k = m
    adjmatrix = np.zeros((N,N),outdtype)
    L0 = 0.5*m0*m
    initialnet = RandomGraph(m0,L0)
    adjmatrix[:m0,:m0] = initialnet
    
    # 2) PERFORM THE PREFERENTIAL ATTACHMENT GROWTH
    # 2.0) Create a list initially containing the hubs ki times
    nodelist = []
    degree = initialnet.sum(axis=1)
    for i in xrange(m0):
        nodelist += [i]*degree[i]
    del degree, initialnet
    
    # 2.1) Incude a new node
    for i in xrange(m0,N):
        counter = 0
        neighbours = []
        while counter < m:
            # 2.2) Choose a node at random
            jid = int(len(nodelist)*numpy.random.rand())
            j = nodelist[jid]
            # Avoid double links
            if adjmatrix[i,j]: continue

            # 2.3) Connect the new node with the randomly chosen node
            adjmatrix[i,j] = 1
            adjmatrix[j,i] = 1
            neighbours.append(j)
            counter += 1

        # 2.4) Update nodelist for new probability in the next iteration
        nodelist += neighbours
        nodelist += m*[i]

    return adjmatrix

def ScaleFreeGraph(N, density, exponent, directed=False):
    """Returns a scale-free graph of given size and exponent.
    
    It follows the method proposed by Goh, Kahng & Kim 'Universal Behaviour
    of Load Distribution in SF networks' PRL 87 (2001). Every node is chosen
    with probability alpha = 1 / (gamma - 1) where gamma is the desired
    exponent.
    
    Parameters
    ----------
    N : integer
        Number of nodes of the final network.
    density : float between 0 and 1.
        Density of links of the graph or the digraph.
    exponent : float
        The exponent (in positive) of the degree-distribution of the resulting
        networks. Recommended values between 2 and 3.
    directed : boolean, optional
        False if a graph is desired, True if digraphs are desired. In case
        of digraphs, both the input and the output degrees follow a scale-
        free distribution but uncorrelated between them.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2 and size NxN.
        The adjacency matrix of the generated scale-free network.
    
    Notes
    -----
    In case of directed networks the input and the output degrees of the
    nodes are correlated, e.g., input hubs are also output hubs.
    """
    # 0) SECURITY CHECKS
    assert density >= 0 and density <= 1, \
        'Density value out of bounds. Values between 0 and 1 accepted.'

    # 1) GENERATE AN EMPTY NETWORK
    adjmatrix = np.zeros((N,N), np.uint8)
    if directed:
        L = round(density*N*(N-1))
    else:
        L = round(0.5*density*N*(N-1))

    # 2) CREATE DEGREE SEQUENCE
    alpha = 1./(exponent - 1.0)
    nodeweights = (np.arange(N) +1)**-alpha
    nodeweights /= nodeweights.sum()    # Probability of a node to be chosen
    nodecumprobability = nodeweights.cumsum()
    del nodeweights

    # 3) CONNECT THE NETWORK
    counter = 1
    while counter <= L:

        # 3.1) Choose two nodes to connect
        xhead = numpy.random.rand()     # A random number between 0 and 1
        xsum = sum(np.sign(nodecumprobability-xhead))
        head = int(0.5*(N-xsum))        

        xtail = numpy.random.rand()
        xsum = sum(np.sign(nodecumprobability-xtail))
        tail = int(0.5*(N-xsum))

        # 3.2) Do not allow self loops and multiple edges
        if head == tail: continue
        if adjmatrix[head,tail]: continue

        # 3.3) If conditions are satisfied, make the link
        adjmatrix[head,tail] = 1
        if not directed:
            adjmatrix[tail,head] = 1
        counter += 1
    
    return adjmatrix

#######################################################################
"""NETWORK REWIRING ALGORITHMS"""
def RewireNetwork(adjmatrix, prewire, directed=False, weighted=False):
    """Returns a network with links rewired with same degrees as net.
    
    The function uses the link switching method to rewire networks while
    conserving the input and output degree of the nodes. See references:
    A.R. Rao & S. Bandyopadhyay, Sankhya, Ser. A 58, 225 (1996),
    J.M. Roberts, Soc. Networks 22, 273 (2000).
    
    Parameters
    ----------
    adjmatrix : ndarray of rank-2
        The adjacency matrix of the network to be rewired.
    prewire : float
        Fraction of links to be rewired. Can be larger than one.
    directed : Boolean, optional
        Specify whether the rewiring gives rise to a directed graph or not.
    weighted : Boolean, optional
        Specify whether the rewiring should only consider binary or weighted
        links. See Usage for more details.
        
    Returns
    -------
    rewmatrix : ndarray of rank-2
    
    Usage
    -----
    1) The parameter 'prewire' is not a probability, but a parameter that
    controls the number of iterations the link swithching will happen. At
    each iteration two links are switched, therefore for a given value, 
    the function will rewire prewire*L links in 1/2*prewire*L iterations,
    where L in the number of links in the network.
    2) If your aim is to create an ensemble of rewired networks for comparison
    of the ensemble properties to those of an empirical network, use some
    value of prewire > 2. We strongly recommend to use prewire = 5 or larger 
    (Zamora-Lopez et al. Structural characterization of networks using the 
    cat cortex as an example. In Lectures in supercomputational
    neuroscience, Springer-Verlag, 2008).
    3) The function creates a list of all the links, requiring more memory.
    The aim of that linklist is to select links with equal probability and
    avoid biases introduced in algorithms that first choose a node at random
    and then one of its links.
    4) If weighted=True, the function will conserve the total input strength
    of the nodes. Conserving both the degrees and the intensities of the nodes
    is not possible, at least one of the four constraints must be broken to
    conserve the other three.
    """
    rewmatrix = np.where(adjmatrix,1,0)
    
    N = len(rewmatrix)
    # 0) GENERATE LIST OF LINKS IF NOT GIVEN.
    if directed:
        linklist = np.array(rewmatrix.nonzero())
    else:
        # Apply nonzero only to the upper triangular part of the matrix
        linklist = np.array(np.triu(rewmatrix).nonzero())

    L = len(linklist[0])
    iterations = int(round(0.5*prewire*L))
    
    count = 0
    while count < iterations:
        # 1) SELECT TWO LINKS:
        linkid1 = int(L*numpy.random.rand())
        linkid2 = int(L*numpy.random.rand())
        # Security check. If the two links are the same, discard the iteration 
        if linkid1 == linkid2: continue
        
        h1 = linklist[0,linkid1]; t1 = linklist[1,linkid1]
        h2 = linklist[0,linkid2]; t2 = linklist[1,linkid2]

        # 2) SECURITY CHECKS
        # Avoid formation of self-loops
        if h1 == t2: continue
        if h2 == t1: continue
        # Avoid formation of double-links
        if rewmatrix[h1,t2]: continue
        if rewmatrix[h2,t1]: continue
        # Avoid trivial exchange of links
        if h1 == h2: continue
        if t1 == t2: continue

        # 3) IF ALL CONDITIONS SUCCESFUL, REWIRE
        # 3.1) Rewire the matrix
        if directed:
            # Put the new links
            if weighted:
                rewmatrix[h1,t2] = rewmatrix[h2,t2]
                rewmatrix[h2,t1] = rewmatrix[h1,t1]
            else:
                rewmatrix[h1,t2] = 1
                rewmatrix[h2,t1] = 1

            # Remove the old links
            rewmatrix[h1,t1] = 0
            rewmatrix[h2,t2] = 0
            
        else:
            # Put the new links
            if weighted:
                rewmatrix[h1,t2] = rewmatrix[h2,t2]; rewmatrix[t2,h1] = rewmatrix[h2,t2]
                rewmatrix[h2,t1] = rewmatrix[h1,t1]; rewmatrix[t1,h2] = rewmatrix[h1,t1]
            else:
                rewmatrix[h1,t2] = 1; rewmatrix[t2,h1] = 1
                rewmatrix[h2,t1] = 1; rewmatrix[t1,h2] = 1

            # Remove the old links
            rewmatrix[h1,t1] = 0; rewmatrix[t1,h1] = 0
            rewmatrix[h2,t2] = 0; rewmatrix[t2,h2] = 0

        # 3.2) update linklist
        linklist[1,linkid1] = t2
        linklist[1,linkid2] = t1

        # 3.3) Count the succesful realization
        count += 1
        
    return rewmatrix

######################################################################
"""MODULAR AND HIERARCHICAL NETWORK MODELS"""
def HMpartition(HMshape):
    """Returns a partition of nodes for a hierarchical/modular random
    network.
    
    Parameters
    ----------
    HMshape : list, tuple or array-like.
        HMshape must be a list of integer numbers indicating the target 
        hierarchical and modular structure. For example, HMshape = [2,4,20]
        will create a network of 2 modules, each module is divided into 4 
        modules of 20 nodes each.
        
    Returns
    -------
    partitions : list
        Contains the partition matrices for each of the nlevel-1 hierarchical
        levels of a network created with given parameter 'HMshape'.
        If HMshape = [2,4,20]
    
    See Also
    --------
    HMRNetwork : Generato of random modular and hierarchical networks.
    """
    N = np.multiply.reduce(HMshape)
    nlevels = len(HMshape)
    
    partitions = []
    for level in xrange(nlevels):
        # 2.1.1) Find the number of blocks per hierarchical level
        if level == 0: continue
        
        nblocks = int(np.multiply.reduce(HMshape[:level]))
        partitionmatrix = np.zeros((N,nblocks), np.uint8)
        
        # 2.2) For each block in the hierarchical level
        for b in xrange(nblocks):
            # 2.2.1) Find the indices for the block
            Nblock = N / nblocks
            iminblock = b*Nblock
            imaxblock = (b+1)*Nblock
            
            partitionmatrix[iminblock:imaxblock,b] = 1
        
        partitions.append(partitionmatrix)

    return partitions

def HMRandomNetwork(HMshape, avklist, directed=False, \
                               outdtype=np.uint8):
    """Returns a random hierarchical/modular network of given shape.
    
    This function generalizes the benchmark hierarchical and modular network
    models introduced in [M.E.J. Newman & M. Girvan, Phys. Rev. E 69, 026113
    (2004)] and in [A. Arenas, A. Diaz-Guilera & C.J. Perez-Vicente, Phys.
    Rev. Lett. 96, 114102 (2006)]. It creates networks of networks that are
    randomly connected at each hierarchical level, i.e., random networks
    that are randomly connected with each other forming larger modules that
    randomly connect with each other. The communities of each level are
    of same size and link density.
    
    Parameters
    ----------
    HMshape : list, tuple or array-like.
        HMshape must be a list of integer numbers indicating the target 
        hierarchical and modular structure. For example, HMshape = [2,4,20]
        will create a network of 2 modules, each module is divided into 4 
        modules of 20 nodes each.
    avklist : list, tuple or array-like.
        The mean degree of the nodes at each hierarchical level. For example,
        avklist = [1,3,15] will generate a network in which nodes will have,
        on average, 15 links with its first hierarchical level community 
        (more internal), 3 links with the rest of communities that form the 
        second hierarchical level (more external), and 1 link with any 
        community of the rest of the network.
    directed : boolean, optional.
        If true, the resulting network will be directed, else, it will be
        undirected.
    outdtype : data type, optional
        Data-type of the resulting adjacency matrix. Default: uint8.
        
    Returns
    -------
    adjmatrix : ndarray of rank-2.
        The adjacency matrix of the generated modular and hierarchical
        network.

    Usage and examples
    ------------------
    The function accepts any desired number of hierarchical levels and
    partitions of the levels into desired number of nodes or communities.
    - HMshape = [4,4,16] will create a network of 4x4x16 = 256 nodes and 
    three hierarchical levels. The network contains 4 communities, each 
    subdivided into 4 communities of 16 nodes.
    - HMshape = [2,5,3,20] will generate a network of 2x5x3x20 = 600 nodes
    and four hierarchical levels. The network is composed of 2 communities, 
    each subdivided into 5 communities, each subdivided into 3 communities 
    of 20 nodes each.
    
    The parameter 'avklist' controls the density of connections between 
    the modules in each of the hierarchical levels.
    - avklist = [1,3,13], for the first example above, means that each node
    is connected with 13 nodes of the first level, to 3 nodes in other
    communities at the second level and makes one link with nodes in the 
    rest of the network which do not belong to the same modules in either
    the first or the second level.
    
    In order to generate the 13-4 or the 15-2 networks in [Arenas et al,
    Phys. Rev. Lett., 114102 (2006)], the function has to be called with
    parameters HMshape = [4,4,16] and avklist = [1,4,13] or avklist = [1,2,15]
    respectively.
    
    Notes
    -----
    The degrees given with parameter 'avklist' are average degrees. Due to
    the random generation process, every node will only have the desired 
    internal and external degrees on the ensemble average, but not
    necessarily in each realization.
    """
    def SeedLinks(adjmatrix, L, partition, directed=False):
        # 0) SSECURITY CHECK
        assert len(partition) > 1, 'Partition needs to have at least 2 communities'
        
        ncoms = len(partition)
        counter = 0
        while counter < L:
            # 1) Choose two communities at random
            com1idx = int(ncoms*numpy.random.rand())
            com2idx = int(ncoms*numpy.random.rand())
            # Make sure the nodes belong to different communities
            if com1idx == com2idx: continue
            
            # 2) Choose one node from each community
            node1 = random.choice(partition[com1idx])
            node2 = random.choice(partition[com2idx])
            # Avoid multiple links
            if adjmatrix[node1,node2]: continue
            
            # Link the nodes and update counter
            if directed:
                adjmatrix[node1,node2] = 1
                counter += 1
            else:
                adjmatrix[node1,node2] = 1
                adjmatrix[node2,node1] = 1
                counter += 1

    # 0) SECURITY CHECKS
    #assert that len(klist) == nlevels
    # 0.2) assert that klist[0] <= n0, klist[1] <= xx, etc.
    
    # 1) PREPARE TO CREATE THE NETWORK
    N = np.multiply.reduce(HMshape)
    nlevels = len(HMshape)
    adjmatrix = np.zeros((N,N), outdtype)
    
    # 2) CREATE THE HM NETWORK BY SEEDING LINKS AT DIFFERENT SCALES
    # 2.1) For each hierarchical level
    for level in xrange(nlevels):
        # 2.1.1) Find the number of blocks per hierarchical level
        if HMshape[level] == 1: continue
        
        if level == 0: nblocks = 1
        else: nblocks = int(np.multiply.reduce(HMshape[:level]))
        
        # 2.2) For each block in the hierarchical level
        for b in xrange(nblocks):
            # 2.2.1) Find the indices for the block
            Nblock = N / nblocks
            iminblock = b*Nblock
            imaxblock = (b+1)*Nblock

            # 2.2.2) Create the partition of nodes into communities within
            # the current block in the current hierarchical level
            ncoms = HMshape[level]
            Ns = Nblock / ncoms
            partition = []
            for p in xrange(ncoms):
                partition.append(range(iminblock + p*Ns,iminblock + (p+1)*Ns))
            
            # 3) SEED THE LINKS AT RANDOM BETWEEN COMMUNITIES
            # NOTE: In the last hierarchical level, every node is one community
            if directed:
                Ls = Nblock * avklist[level]
            else:
                Ls = 0.5 * Nblock * avklist[level]
            
            SeedLinks(adjmatrix, Ls, partition, directed)
    
    return adjmatrix
