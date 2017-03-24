#!/usr/bin/env python
# Setup Script
# Author: Leon Foks
# March 10 2017
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
      package_data={'gatdaem1d':['gatdaem1d.so']},
      scripts=[],
      version=1.0,
      description='Time domain electromagnetic forward modeller',
      long_description='',
      classifiers=[
        'Development Status :: 0 - Alpha',
        'License :: OSI Approved :: GPL?',
        'Programming Language :: Python :: 3.5',
        'Topic :: Electromagnetic :: forward modelling :: geophysics',
      ],
      author='Ross Brodie',
      author_email='ross.c.brodie@ga.gov.au',
      install_requires=[
          'numpy>=1.11',
      ],
      url='https://github.com/GeoscienceAustralia/ga-aem')
