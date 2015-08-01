TODO list for GAlib
===================


High priority...
----------------

#. Identify and modify all functions whose unweighted version could be improved using Boolean operators: ``ReciprocalDegree()``, ``RichClub()``, ``MatchingIndex()``, ``ConnectedComponents()``, ``K_Core()``, ``K_Shells()``, ``Modularity()``, etc.
#. Identify the slowest functions which could significantly be accelarated using Numba package and create Numba-based duplicates.
#. Write a function for the Dijkstra algorithm.
#. I/O support to more graph formats: graphML (.xml), DOT (.dot), etc. 
#. Write tests for all methods
#. Replace dangerous in-place modification functions by pure functions (returning a new object)
#. *Suggest your own.*

Would be nice to ...
--------------------

* Finish documentation (use Sphinx for that).
* Check compatibility with scipy.sparse module.
* *Suggest your own.*
