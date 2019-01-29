## TODO list for GAlib


### Priority...


1. Identify and modify all functions whose unweighted version could be improved using Boolean operators: ``ReciprocalDegree()``, ``RichClub()``, ``MatchingIndex()``, ``ConnectedComponents()``, ``K_Core()``, ``K_Shells()``, ``Modularity()``, etc.
2. Identify further functions which could be accelarated using Numba package, and write the Numba-based duplicates.
3. Write a function for the Dijkstra algorithm.
4. Support for weighted network metrics (ACHTUNG!! I won't accept just any algorithm for weighted networks, since many do not make sense.
4. I/O support to more graph formats: igraph, graphML (.xml), DOT (.dot), etc. 
5. Write tests for all methods
6. *Please, suggest your own…*

### Would be nice to ...

* Finish documentation (use Sphinx for that).
* Check compatibility with scipy.sparse module.
* *Please, suggest your own…*
