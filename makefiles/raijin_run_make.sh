#!/bin/bash

#Set the directory path for dependencies
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export tntdir='../submodules/tnt'

#GNU compiler on raijin.nci.org.au
module load ga-aem/dev-gnu
export cxx=g++
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always'
export exedir='../bin/raijin/gnu'

#Intel compiler on raijin.nci.org.au
#module load ga-aem/dev-intel
#export cxx=icpc
#export mpicxx=mpiCC
#export cxxflags='-std=c++11 -O3 -Wall -diag-disable remark'
#export exedir='../bin/raijin/intel'


echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme
echo ---------------------------------------


#Compiled as shared libs
make -f gatdaem1d_python.make $1
make -f gatdaem1d_matlab.make $1

#Compile without MPI
make -f ctlinedata2sgrid.make $1
make -f ctlinedata2slicegrids.make $1
make -f example_forward_model.make $1
make -f gaforwardmodeltdem.make $1

#Compile with MPI
make -f galeisbstdem.make $1
make -f garjmcmctdem.make $1
make -f galeiallatonce.make $1
make -f galeisbsfdem.make $1



