## TODO list for GAlib


### Priorities...

1. Drop support for Python 2.
2. Clean-up the library files. Remove unnecessary comments, copyright duplicates, etc.
3. Update to the newer packaging and PyPI release standards.
4. Integrate PathLims into pyGAlib.
5. Bring weighted network generation and randomization from SiReNetA.
6. Add the generation of random graphs with specified degree(-degree) correlations. 


### This and that...

1. Identify further functions which could be accelarated using Numba package, and write the Numba-based duplicates.
2. Write a function for the Dijkstra algorithm.
3. Wrtie a function to calculate the Katz centrality.
3. Identify and modify all functions whose unweighted version could be improved using Boolean operators: ``ReciprocalDegree()``, ``RichClub()``, ``MatchingIndex()``, ``ConnectedComponents()``, ``K_Core()``, ``K_Shells()``, ``Modularity()``, etc.
4. Support for weighted network metrics (ACHTUNG!! I won't accept just any algorithm for weighted networks, since many do not make sense.)
5. I/O support to more graph formats: igraph, graphML (.xml), DOT (.dot), etc. 
6. Write tests for all functions.
7. *Suggest your own…*


### It would be nice to ...

1. Add functions for community detection (please help!).
2. Add functionalities for graph visualization.
3. Finish documentation (use Sphinx for that).
4. Check compatibility with scipy.sparse module.
5. *Please, suggest your own…*
