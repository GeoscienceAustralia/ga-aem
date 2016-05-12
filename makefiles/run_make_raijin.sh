#!/bin/sh

#Script to load compiler modules and dependent software 
#on raijin.nci.org.au

module load gcc/5.2.0
module load openmpi/1.6.3
module load fftw3/3.3.3

#make -f gaforwardmodeltdem.make
#make -f galeisbstdem.make
#make -f garjmcmctdem.make

#make -f gatdaem1d_matlab.make allclean
make -f gatdaem1d_python.make allclean

