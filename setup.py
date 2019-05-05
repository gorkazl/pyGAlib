"""setup.py"""

from setuptools import setup, find_packages


with open("requirements.txt") as reqs_file:
    REQS = [line.rstrip() for line in reqs_file.readlines() if line[0] not in ['\n', '-', '#']]

config = {
    'name': 'galib',
    'description': ('A library for graph analysis in Python / NumPy.'),
    # Check these exist in PyPI
    'keywords': 'graph theory, complex networks, network analysis',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    'author': 'Gorka Zamora-Lopez',
    'author_email': 'galib@Zamora-Lopez.xyz',
    'url': 'https://github.com/gorkazl/pyGAlib',
    'version': '1.0.2',
    'license': 'Apache License 2.0',
    'install_requires': REQS,
    'packages': find_packages(exclude=['doc', '*tests*']),
    'scripts': [],
    'include_package_data': True
    }

setup(**config)
