#!/bin/sh

#Script to load compiler modules and dependent software 
#on raijin.nci.org.au

module load gcc/5.2.0
module load openmpi/1.6.3
module load fftw3/3.3.3
module load python3/3.3.0

make -f gaforwardmodeltdem.make allclean
#make -f galeisbstdem.make allclean
#make -f garjmcmctdem.make allclean

#make -f gatdaem1d_matlab.make allclean
#make -f gatdaem1d_python.make allclean

