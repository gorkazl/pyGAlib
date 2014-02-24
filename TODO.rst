TODO list for GAlib
===================


High priority...
----------------

1. Identify and modify all functions whose unweighted version could be improved using Boolean operators:

    ``ReciprocalDegree()``, ``RichClub()``, ``MatchingIndex()``, ``ConnectedComponents()``, ``K_Core()``, ``K_Shells()``, ``Modularity()``.

2. Identify the slowest functions which could significantly be accelarated using Numba package and create Numba-based duplicates.
3. Write a function for the Dijkstra algorithm.
4. I/O support to more graph formats: graphML (.xml), DOT (.dot), etc. 
5. Suggest your own.

Would be nice to ...
--------------------

* Finish documentation (use Sphinx for that).
* Check compatibility with scipy.sparse module.
* Suggest you own.
