#!/bin/sh

#Script to load compiler modules and dependent software 

module load fftw

#GNU compiler on raijin.nci.org.au
export cxx=g++
export cc=gcc
export cxxflags='-std=c++11 -O3 -Wall'

make -f gatdaem1d_python.make allclean

