#!/bin/sh

module load cmake/3.16.2
module load fftw3/3.3.8
module load netcdf/4.8.0
module load gdal/3.0.2
module load openmpi/3.0.4
module load petsc/3.12.2

module load intel-compiler
#module load gcc/11.1.0

export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
echo $PKG_CONFIG_PATH

export BUILD_DIR=$PWD/build-gnu
export INSTALL_DIR=$PWD/install-gnu

module list 

# build is a temporary directory for building (compiling and linking)
mkdir $BUILD_DIR
cd $BUILD_DIR

#eg for the GNU compilers
#cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
#eg for the Intel compilers
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release  ..


# Build all (cmake --build .)
cmake --build . --target all
#cmake --build . --target galeisbstdem.exe
#cmake --build . --target galeiallatonce.exe

# Build all
#cmake --build . --target clean

# Build only specific items
#cmake --build . --target cpp-utils-static
#cmake --build . --target galeisbstdem.exe

# Install executables and libraries into $INSTALL_DIR
cmake --install . --prefix $INSTALL_DIR

