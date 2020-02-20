#!/bin/sh

#Script to load compiler modules and dependent software 

module load openmpi/1.6.3
module load fftw3/3.3.3

#GNU compiler on raijin.nci.org.au
module load gcc/5.2.0
export cxx=g++
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall'
export exedir='../bin/vdi/gnu'

export cpputilssrc='../submodules/cpp-utils/src'
export tntdir='../submodules/tnt'
export srcdir='../src'

#Intel compiler on raijin.nci.org.au
#module load intel-cc/14.1.106
#export cxx=icpc
#export mpicxx=mpiCC
#export cxxflags='-std=c++11 -O3 -Wall -diag-disable remark'
#export exedir='../bin/raijin/intel'

#make -f gaforwardmodeltdem.make allclean
#make -f galeisbstdem.make allclean
make -f garjmcmctdem.make allclean
#make -f gatdaem1d_python.make allclean
#make -f gatdaem1d_matlab.make allclean

