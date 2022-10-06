#!/bin/sh

# BUILD_DIR is a temporary directory for building (compiling and linking)
export BUILD_DIR=$PWD/build-ubuntu
# INSTALL_DIR is the directory for installing the build package
export INSTALL_DIR=$PWD/install-ubuntu

mkdir $BUILD_DIR
cd $BUILD_DIR

#Example for the GNU compilers
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --target all
cmake --install . --prefix $INSTALL_DIR

#cmake --build . --target galeisbstdem
#cmake --build . --target garjmcmctdem
#cmake --build . --target galeisbsfdem
#cmake --build . --target galeiallatonce

#cmake --build . --target gaforwardmodeltdem
#cmake --build . --target example_forward_model
#cmake --build . --target example_forward_model_c
#cmake --build . --target gatdaem1d-static
#cmake --build . --target gatdaem1d-shared

#cmake --build . --target ctlinedata2sgrid
#cmake --build . --target ctlinedata2slicegrids
#cmake --build . --target removelog10conductivityfromsgrid
## Windows only ctlinedata2georefimage
## Windows only ctlinedata2curtainimage


