# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'Gorka Zamora-López (Ph.D.)',
    'url': 'https://github.com/gorkazl/pyGAlib',
    'download_url': 'https://github.com/gorkazl/pyGAlib/archive/master.zip',
    'author_email': 'Gorka.zamora@ymail.com',
    'version': '0.1',
    'install_requires': ['numpy'],
    'packages': ['pyGAlib'],
    'scripts': [],
    'name': 'pyGAlib'
}

setup(**config)
