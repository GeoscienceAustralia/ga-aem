#!/bin/sh

#Script to load compiler modules and dependent software 
#on raijin.nci.org.au

module load gcc/5.2.0
module load fftw3/3.3.3

make -f gatdaem1d_matlab.make allclean

