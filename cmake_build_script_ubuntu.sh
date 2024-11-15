#!/bin/sh

# INSTALL_DIR is the directory for installing the build package
export INSTALL_DIR=$PWD/install-ubuntu

# BUILD_DIR is a temporary directory for building (compiling and linking)
export BUILD_DIR=$PWD/build-ubuntu

mkdir $BUILD_DIR
cd $BUILD_DIR

# Switches for turning off certain dependencies if they are not wanted or available
# 	-DWITH_FFTW=OFF
# 	-DWITH_MPI=OFF
# 	-DWITH_NETCDF=OFF
# 	-DWITH_GDAL=OFF
# 	-DWITH_PETSC=OFF
# cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF -DWITH_NETCDF=OFF -DWITH_GDAL=OFF -DWITH_PETSC=OFF ..

#Example for the GNU compilers
cmake -Wno-dev -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
cmake --build . 
cmake --install . --prefix $INSTALL_DIR

# Or alternatively ...

# Build only particular targets
#cmake --build . --target galeisbstdem
#cmake --build . --target galeisbstdem-nompi
#cmake --build . --target garjmcmctdem
#cmake --build . --target galeiallatonce

#cmake --build . --target gaforwardmodeltdem
#cmake --build . --target example_forward_model
#cmake --build . --target example_forward_model_c

#cmake --build . --target gatdaem1d-static
#cmake --build . --target gatdaem1d-shared
#cmake --build . --target matlab-bindings --config=Release
#cmake --build . --target python-bindings --config=Release

#cmake --build . --target ctlinedata2sgrid
#cmake --build . --target ctlinedata2slicegrids
#cmake --build . --target removelog10conductivityfromsgrid
## Windows only ctlinedata2georefimage
## Windows only ctlinedata2curtainimage



