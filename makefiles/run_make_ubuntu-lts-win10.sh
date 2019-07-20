#!/bin/bash

#Set the directory path for dependencies
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export tntdir='../submodules/tnt'

#GNU compiler in ubuntu 18.04 LTS app on Windows 10 
export cxx=g++
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wno-unused-result -Wno-format-security -fdiagnostics-color=always'
export exedir='../bin/ubuntu/gnu'
export FFTW_DIR='/usr/lib/x86_64-linux-gnu'
export PETSC_DIR='/usr/lib/petscdir/3.7.7/x86_64-linux-gnu-real'

echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme
echo ---------------------------------------

#Compiled as shared libs
#make -f gatdaem1d_python.make $1
#make -f gatdaem1d_matlab.make $1

#Compile without MPI
make -f ctlinedata2sgrid.make $1
#make -f ctlinedata2slicegrids.make $1
make -f example_forward_model.make $1
make -f gaforwardmodeltdem.make $1

#Compile with MPI
make -f galeisbstdem.make $1
make -f garjmcmctdem.make $1
make -f galeiallatonce.make $1
make -f galeisbsfdem.make $1


