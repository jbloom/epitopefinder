"""Setup script for ``epitopefinder``.

This script uses ``distutils``, the standard python mechanism for installing
packages. To build, test, and install the package, use the following
commands::

    python setup.py build
    python setup.py test
    python setup.py install

If the user does not have permissions to write to the install directory,
the last command may need to be replaced by::

    sudo python setup.py install

In order for plotting to be enabled, ``pylab`` and ``matplotlib`` must be
installed and available. If they are not available, this script prints a
warning indicating that fact.

Written by Jesse Bloom.
"""


import sys
import os
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command


# main setup command
setup(
    name = 'epitopefinder', 
    version = '0.1', 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fhcrc.org', 
    url = 'https://github.com/jbloom/epitopefinder', 
    description = 'Finds immune epitopes and processes them',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Mac OS X.',
    packages = ['epitopefinder'],
    package_dir = {'epitopefinder':'src'},
    scripts = [
            'scripts/epitopefinder_getepitopes.py',
            'scripts/epitopefinder_plotlineardensity.py',
            'scripts/epitopefinder_plotcorrelation.py',
            'scripts/epitopefinder_plotdistributioncomparison.py',
            'scripts/epitopefinder_selectsites.py',
            ],
)
