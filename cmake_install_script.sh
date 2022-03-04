#!/bin/sh

# build is a temporary directory for building (compiling and linking)
mkdir build
cd build

#eg for the GNU compilers
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
#eg for the Intel compilers
#cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release  ..

# Build all (cmake --build .)
cmake --build . --target all

# Build all
#cmake --build . --target clean

# Build only specific items
#cmake --build . --target cpp-utils-static
#cmake --build . --target galeisbstdem.exe

# Install executables and libraries into my_install_dir
#cmake --install . --prefix ../my_install_dir

