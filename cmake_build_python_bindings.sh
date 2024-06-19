#!/bin/sh

#module load gcc/13.2.0
#module load cmake/3.21.4
#module load fftw3/3.3.8
#module list 

# Add the PETSc pkg-config path to the pkg-config search path
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# pkg-config for PETSc seems non-standard on gadi so add the ompi/compiler specific library directory that contains the library file libpetsc.so
export PETSC_LIBRARY_DIR=$PETSC_DIR/lib/ompi3/GNU

# BUILD_DIR is a temporary directory for building (compiling and linking)
export BUILD_DIR=$PWD/build-gnu
# INSTALL_DIR is the directory for installing the build package
export INSTALL_DIR=$PWD/install-gnu

mkdir $BUILD_DIR
cd $BUILD_DIR

cmake -Wno-dev -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF -DWITH_NETCDF=OFF -DWITH_GDAL=OFF -DWITH_PETSC=OFF -DCMAKE_PREFIX_PATH=$FFTW_ROOT .. 
cmake --build . --target python-bindings
cmake --install . --prefix python/gatdaem1d
