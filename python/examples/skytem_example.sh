#!/bin/bash

#Make sure you are using Python3.x (Python2.x will not work)
#Make sure the FFTW libraries are in your search path

#I have found that the gcc used to compile the shared library also needs to be in the path

#Make sure the directory --/--/ga-aem/python is in your PYTHONPATH environment varible
#This will fail: "export PYTHONPATH='blah/blah/ga-aem/python/gatdaem1d"
#This will work: "export PYTHONPATH='blah/blah/ga-aem/python"

#module load gcc/5.2.0
#module load fftw3/3.3.3
#module load python3/3.3.0
#module load python3/3.3.0-matplotlib


export PYTHONPATH='../../python'
python3 skytem_example.py

