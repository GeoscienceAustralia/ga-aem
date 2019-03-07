#!/bin/sh

#Script to load compiler modules and dependent software 
module load gcc/5.2.0
module load openmpi/1.6.3
module load fftw3/3.3.3   #required for forward modelling and inversion
module load petsc/3.4.3   #required for galeiallatonce only
module load gdal/1.11.1   #required for ctlinedata2slicegrids only
module load intel-mkl/14.1.106 #required for galeiallatonce only

#Set the directory path for dependencies
export cpputilssrc='../submodules/cpp-utils/src'
export tntdir='../submodules/tnt'
export srcdir='../src'

#GNU compiler on virtual machine
export cxx=g++
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always'
export exedir='../bin/vm/gcc'

make -f example_forward_model.make $1
make -f gaforwardmodeltdem.make $1
make -f gatdaem1d_python.make $1
make -f gatdaem1d_matlab.make $1

make -f galeisbstdem.make $1
make -f garjmcmctdem.make $1
make -f galeiallatonce.make $1

make -f ctlinedata2sgrid.make $1
make -f ctlinedata2slicegrids.make $1

make -f galeisbsfdem.make $1

