#!/bin/bash

#module load fftw3/3.3.3

#export PYTHONPATH=/short/public/rcb547/apps/gaaem-1.0/python
export PYTHONPATH=/short/public/rcb547/apps/gaaem-1.0/python/gatdaem1d

#export PATH=$PATH:/short/public/rcb547/apps/gaaem-1.0/python/bin/raijin/gnu
#export PATH=$PATH:/short/public/rcb547/apps/gaaem-1.0/python/gatdaem1d
#echo $PATH

python skytem_example.py

