## TODO list for GAlib


### Priorities...

1. ~~Drop support for Python 2.~~
2. ~~Clean-up the library files. Remove unnecessary comments, copyright duplicates, etc.~~
3. ~~f" … " string formatting~~. **Was not needed, actually**.
4. ~~Replace runtime "prints" by proper warning and error detection.~~
5. Update to the newer packaging and PyPI release standards.
    1. We will use Hatch, at least for now.
    2. Move Matploblib into optional dependencies. Only used to run the examples.
    3. Prepare GAlib for conda-sourceforge.
    4. Add in README.md instructions for conda users (either `conda install galib` or, first install dependencies via conda and then `python -m pip install -e . --no-deps`. Or, release also a *yml* file with preinstallation of the dependencies.
6. ~~What should I do with the Examples/ and Branding/ folders when packaging?~~ Should I integrate them to the wheel, or should I let users to download the examples separately (manually) from the GitHub page? **NO, do not include them into the wheel. Just leave them in the root of the repo, for independent download.**
7. Clean, revise and test scripts in the Examples/ folder.
8. Integrate PathLims into pyGAlib (?).
9. Bring weighted network generation and randomization from SiReNetA.
10. Add the generation of random graphs with specified degree(-degree) correlations. 
11. Update README.md file:
    1. ~~Add version and status indicators at the top~~. 
    2. ~~Revise installation instructions.~~



### This and that...

1. Identify further functions which could be accelarated using Numba package, and write the Numba-based duplicates.
2. Write a function for the Dijkstra algorithm.
3. Write a function to calculate the Katz centrality.
4. Add further data I/O conversions to other graph formats: igraph, graphML (.xml), DOT (.dot), etc. 
5. Identify and modify all functions whose unweighted version could be improved using Boolean operators: maybe ``ReciprocalDegree()``, ``RichClub()``, ``MatchingIndex()``, ``ConnectedComponents()``, ``K_Core()``, ``K_Shells()``, ``Modularity()``, etc.
6. Support for weighted network metrics (ACHTUNG!! I won't accept just any algorithm for weighted networks, since many do not make sense.)
7. Write tests for all functions.
8. Finish and publish documentation.
9. *Suggest your own…*


### It would be nice to ...

1. Use some *linting* software to double check code style.
2. Add functions for community detection (*please help!*).
3. Add functionalities for graph visualization (*please help!*).
5. Check compatibility with scipy.sparse module.
6. *Please, suggest your own…*
