"""setup.py"""

from setuptools import setup, find_packages


with open("requirements.txt") as reqs_file:
    REQS = [line.rstrip() for line in reqs_file.readlines() if line[0] not in ['\n', '-', '#']]

setup(
    name = 'galib',
    description = 'A library for graph analysis in Python / NumPy.',
    version = '1.0.4',
    url = 'https://github.com/gorkazl/pyGAlib',
    license = 'Apache License 2.0',

    author = 'Gorka Zamora-Lopez',
    author_email = 'galib@Zamora-Lopez.xyz',

    install_requires = REQS,
    packages = find_packages(exclude=['doc', '*tests*']),
    scripts = [],
    include_package_data = True,

    keywords = 'graph theory, complex networks, network analysis',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
