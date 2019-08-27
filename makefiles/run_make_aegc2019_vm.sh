#!/bin/sh

#Script to load compiler modules and dependent software 

#Must enable the devtoolset-4 first
#scl enable devtoolset-4 bash


#Set the directory path for dependencies
export cpputilssrc='../submodules/cpp-utils/src'
export tntdir='../submodules/tnt'
export srcdir='../src'
export FFTW_DIR='/usr/lib64'

#GNU compiler on virtual machine
export cxx=g++
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always -I/usr/include/gdal'
export exedir='../bin'

#make -f example_forward_model.make $1
#make -f gaforwardmodeltdem.make $1
#make -f gatdaem1d_python.make $1
#make -f gatdaem1d_matlab.make $1
#make -f galeisbstdem.make $1
make -f garjmcmctdem.make $1
#make -f galeiallatonce.make $1

#make -f ctlinedata2sgrid.make $1
#make -f ctlinedata2slicegrids.make $1

#make #-f galeisbsfdem.make $1



