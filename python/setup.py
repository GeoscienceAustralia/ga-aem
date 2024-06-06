#!/usr/bin/env python
# Setup Script
# Author: Leon Foks
# March 10 2017
# Updated: Ross C Brodie, March 20 2017
# Updated: Ross C Brodie, March 21 2024 to Version 2.0
# Updated: Leon Foks, June 5 2024. Added .dylib 
import sys
import os
from os.path import join

# Test Python's version
major, minor = sys.version_info[0:2]
if (major, minor) < (3, 5):
    sys.stderr.write('\nPython 3.5 or later is needed to use this package\n')
    sys.exit(1)

try:
    from setuptools import setup
except ImportError:
    pass

setup(name='gatdaem1d',
      packages=['gatdaem1d'],
      package_dir={'gatdaem1d':'gatdaem1d'},
      package_data={'gatdaem1d':['gatdaem1d.so', 'gatdaem1d.dll', 'gatdaem1d.dylib']},
      scripts=[],
      version=2.0,
      description='Time-domain airborne electromagnetic forward modelling.',
      long_description='Time-domain airborne electromagnetic forward modelling. Python interface to C++ library for forward model and derivative calculations for airborne electromagnetic (AEM) systems used in geophysics.',
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Electromagnetic :: Airborne :: Forward modelling :: Geophysics',
      ],
      author='Ross C Brodie, Geoscience Australia and Leon Foks',
      author_email='ross.c.brodie at ga.gov.au',
      install_requires=[
          'numpy>=1.11',
      ],
      url='https://github.com/GeoscienceAustralia/ga-aem',
      license='GNU General Public License v2 (GPLv2)')
