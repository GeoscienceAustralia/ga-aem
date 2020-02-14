#!/bin/sh

#export NCI_GXX_ABI_WARNING=1

#Set the directory path for dependencies
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'

#GNU compiler on raijin.nci.org.au
export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always -D_GLIBCXX_USE_CXX11_ABI=1 -Wno-unknown-pragmas'
module load gcc/5.2.0
module load fftw3/3.3.7-gcc #Use the same module as Julia just to be sure
module list

rm gatdaem1d_julia.so
#Compile as shared lib
g++ -fPIC -shared -O3 $cxxflags \
-I$srcdir \
-I$cpputilssrc \
-L$FFTW_DIR -lfftw3 \
$cpputilssrc/general_utils.cpp \
$cpputilssrc/file_utils.cpp \
gatdaem1d_julia.cpp \
-o gatdaem1d_julia.so


