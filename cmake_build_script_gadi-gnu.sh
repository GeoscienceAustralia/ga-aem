#!/bin/sh

module load gcc/13.2.0
module load cmake/3.21.4
module load fftw3/3.3.8
module load netcdf/4.8.0
module load gdal/3.0.2
module load openmpi/4.1.4
module load petsc/3.17.4
module list 

# Add the PETSc pkg-config path to the pkg-config search path
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# pkg-config for PETSc seems non-standard on gadi so add the ompi/compiler specific library directory 
export PETSC_EXTRA_LIB_DIR=$PETSC_DIR/lib/ompi3/Intel
#export PETSC_EXTRA_LIB_DIR=$PETSC_DIR/lib/ompi3/GNU

# BUILD_DIR is a temporary directory for building (compiling and linking)
export BUILD_DIR=$PWD/build-gnu
# INSTALL_DIR is the directory for installing the build package
export INSTALL_DIR=$PWD/install-gnu

mkdir $BUILD_DIR
cd $BUILD_DIR

cmake -Wno-dev -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
cmake --build . --target all
cmake --install . --prefix $INSTALL_DIR

#cmake --build . --target galeisbstdem
#cmake --build . --target galeisbstdem-nompi
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

