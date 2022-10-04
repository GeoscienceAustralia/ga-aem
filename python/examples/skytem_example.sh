#!/bin/bash

#Make sure you are using Python3.x (Python2.x will not work)
#Make sure the FFTW libraries are in your search path

#Make sure the directory <ga-aem-install-dir>/python is in your PYTHONPATH environment variable
# -This will fail: "export PYTHONPATH='blah/blah/ga-aem/python/gatdaem1d"
# -This will work: "export PYTHONPATH='blah/blah/ga-aem/python"

#I have found that the gcc used to compile the shared library also needs to be in the path

# You may need to load some modules on gadi.nci.org.au
#module load gcc/11.1.0
#module load fftw3/3.3.8
#module load python3/3.10.4

export PYTHONPATH='../../python'
python3 skytem_example.py

