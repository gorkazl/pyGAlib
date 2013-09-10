GRAPH ANALYSIS LIBRARY (v0.0 beta)
(A preliminary manual)
by Dr. Gorka Zamora-López



WHAT IS GAlib?

GAlib is a library for the analysis of graphs and complex networks in Python/NumPy. It is intended for researchers studying complex networks and those performing data analysis with a reasonable knowledge of graph theory and Python/NumPy programming. The library is very easy to modify and extend in order to adapt it for the users' personal needs. As a network scientist and developer of data analysis methods myself, flexibility was the original motivation to write GAlib. It does not include visualization functionalities but it can be easily integrated in larger projects.

GAlib does not define a graph object, it represents networks by their adjacency matrices as rank-2 ndarrays. This choice limits the size of the networks GAlib can handle but it allows to exploit NumPy to boost performance far beyond code written in pure Python.

The library contains 3 modules. See the Reference Guide or the interactive documentation for a detailed list of functions within each module.
	- gatools.py : miscelaneous helper functions, e.g. data-type transformations.
	- galib.py : basic graph descriptors, e.g. degrees and clustering.
	- gamodels.py : Generation of synthetic networks and randomization.

In addition to the library files GAlib comes with several example scripts to illustrate its use (see .../GAlib/Examples/ folder). For network scientists who are not familiar with the Python / NumPy environment, GAlib includes also a few helper scripts to generate and save into files ensembles of random or scale-free networks, as-well-as ensembles of rewired networks conserving the degrees. Users only need to modify a few selected fields and let GAlib do the rest. See .../GAlib/HelperScripts/ folder.


INSTALLATION

As a library, GAlib does not require any installation but the modules need to be placed in the python path such that they can be imported. The only requirement is to have NumPy/SciPy installed together with Python. If NumPy/Scipy is not yet installed in the computer, visit the following webpage for further information (http://www.scipy.org).

For beginners to Python/Numpy/Scipy it is highly recommended to install a Python distribution that already includes NumPy, SciPy and other utilities for numerical and scientific tools. These are easy to install and stable distributions. For example Canopy (https://www.enthought.com/products/canopy/) or Anaconda (https://store.continuum.io/cshop/anaconda/)

If you don't aim to modify the code in GAlib the simplest way to install it is to copy the modules into the site-packages folder. This is a special folder that is installed with Python to place third party libraries and packages. However, its exact location depends both on the operating system and the Python version(s) or distribution(s) currently installed. 

	1) Find the current python version installed.
	This might be more complicated and annoying than it sounds because there might be several versions of Python installed in the computer. For example, Mac OS X and Linux come with a pre-installed version of Python. When re-installing Python through a distribution the older version is not removed and each version has its own "site-packages" folder. To find the current Python version open a terminal and type "python". An interactive shell will open (type exit() to leave the interactive shell). The current Python version is displayed in the header. If a Python distribution has been installed that includes NumPy and SciPy, very likely iPython has been installed too. iPython is an advanced interactive shell for Python. If this is the case, open another terminal and type ipython to open another interactive shell. In the header the version of Python used by iPython is displayed. Make sure that the two version displayed in the headers are the same. If the two versions are different you will have to copy GAlib into the site-packages folder of both Python versions to warranty that it will always work independently of how to run the scripts.

	2) Locate the "site-packages" folder.
	Once the current version of installed Python installed is known, search the filesystem for folders named "site-packages". Usually the whole path looks like "…/lib/pythonX.Y/site-packages" where X.Y is the Python version number. If several Python versions are installed in the computer, locate the "site-packages" folders for those versions you will be using.

	3) Copy the modules of GAlib into the "site-packages" folder.

	4) Check GAlib is correctly installed.
	Close any open interactive shells and open them again by typing python or ipython in the terminal. Import the modules by typing import galib, import gatools and import gamodels. If no error is returned, the library is correctly installed. Otherwise a warning will be raised that the modules are not found.

Now, if you intend to modify the library or to include your own functions to it, the modules in the "site-packages" are usually not possible to modify due to permission restrictions. Alternatively, you can place the modules of GAlib in any folder you wish, e.g. "/usr/myusername/GraphAnalysisLibrary". To allow Python to find GAlib a text file can be created in the site-packages folder that contains additional folders to be included in the path. Follow steps 1) and 2) previously and do the following:

	3) Include a .pth file in the site-packages folder
Open a new text document in an editor and type in the files the path you want to include, e.g. the file will have a single line with the text /usr/myusername/GraphAnalisysLibrary. Save the text document in the site-packages folder and give it the extension .pth, for example: extrapythonpaths.pth. You can include as many paths as you want in a .pth file, every path being in one independent line.



FUTURE PLANS

At this starting stage GAlib includes basic graph analysis tools. In future releases further functionalities, graph measures and network model generators will be included:
	- measures for weighted and/or directed networks,
	- classical graph models,
	- functions to compute the roles of nodes in networks with modular organization,
	- improved data conversions for graph formats used in other packages (Pajek, graph-tools, NetworkX, Gephi, etc.)
	- support for sparse matrices to increase the size of networks handled.

Any third party collaboration to improve and extend GAlib is highly welcome.


1.4) HOW TO FIND DOCUMENTATION
While working in an interactive session type:

>>> help(modulename)

where modulename can be galib, gamodels or gatools, to see a list of functions included in the module. To find further information of each individual function: aim, parameters and returned values, type:

>>> help(functionname)

In IPython the help command is replaced by a question mark after the module or function name:

>>> modulename?
>>> functionname?



WHERE IT ALL COMES FROM?

GAlib is an accidental by-product of my scientific research. In 2004 I started my Ph.D. in the group of Prof. Jürgen Kurths at the University of Potsdam (Germany) to make research in brain connectivity. My goal was to uncover the topological organization of long-range connections in the cat and the macaque brains. Back in those times complex networks was a exploding field and the computational tools were scarce. There were some very good programs, e.g. Pajek, but such pre-packaged programs were too rigid environments for a field constantly in change. I couldn't just wait until someone would include new graph measures and models to those programs and I just started coding myself every graph measure that I needed.  Hence, the need for flexibility in a research environment was the main reason for GAlib to happen, and remains its main scope. GAlib is not intended for unexperienced users that need a "blackbox" to put their data in. I enjoyed that process very much and I ended making my own contributions to graph theory, particularly on the meaning and quantification of significance of graph measures.

Over the years I was often tempted to publicly release the library and so contribute to the open-source ecosystem of Python and NumPy. The evolution of Numeric and Numarray into NumPy delayed my intentions, and it was for good. The extra experience I have gained in the mean time allowed me to improve every single function from code that "just works" to code that I hope to be very much optimized. I discovered that optimizing Python code is very non-trivial. Far from the mandates given by the Zen of Python, the combination of Python and NumPy allows to very often do the same job in several different manners, each having its own pros and cons. in Appendix II I summarize some of those "tricks".


