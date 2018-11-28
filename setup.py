'''setup.py'''

from setuptools import setup, find_packages


with open("requirements.txt") as reqs_file:
    REQS = [line.rstrip() for line in reqs_file.readlines() if line[0] not in ['\n', '-', '#']]

config = {
    'name': 'galib',
    'description': ('Write me.'),
    'keywords': 'hbp, human brain project, collaboratory, library, science',
    'classifiers': [
        'Development Status :: 0 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    'author': 'HBP Infrastructure Developers',
    'author_email': 'platform@humanbrainproject.eu',
    'url': 'https://github.com/HumanBrainProject/hbp-service-client',
    'version': '0.0.1',
    'license': 'Apache License 2.0',
    'install_requires': REQS,
    'packages': find_packages(exclude=['doc', '*tests*']),
    'scripts': [],
    'include_package_data': True
    }

setup(**config)
