# GA-AEM Source Code Repository
[![Build Status](https://travis-ci.com/GeoscienceAustralia/ga-aem.svg?branch=develop)](https://travis-ci.com/GeoscienceAustralia/ga-aem)

## Geoscience Australia Airborne Electromagnetics Programs

- Authors
	- Ross C Brodie, Geoscience Australia (ross.c.brodie at ga.gov.au)
	- Richard Taylor, Geoscience Australia

- Language: Mostly C++, some matlab, some python

## Cloning

When initially cloning the repository in git it is convenient to use the "--recursive" option so that all the submodules and their respective submodules are initialised and populated with code.

```bash
> git clone --recursive git@github.com:rcb547/ga-aem-csiro.git
```

## Creation Notes
- This repo was duplicated from https://github.com/GeoscienceAustralia/ga-aem.git on 20210706 
- It was as made private for RB's CSIRO Secondment work but allows eventual pushing back to the original repo.
- Followed instructions from the link below
	- https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/duplicating-a-repository

## Releases
### Release-20160606
	- Added Python 3.x interface for simple forward modelling and derivatives only.
	- Added Matlab interface for simple forward modelling and derivatives only.
	- Changed how the PPM normalisation is carried out. Now PPM normalisation is by directional-component-wise with respect to the maximum primary dB/dt or B-field at the receiver for a reference system geometry (which is usually estimated on a per flight or per survey basis). Previously PPM normalisation was with respect to the system geometry for the forward model being run.
	- Added GEOTEM (1996 ppm system) and SPECTREM (ppm system) examples.
	- Fixed a bug in the thickness derivative of the second bottom layer. This may have effected few-layer inversions, but not multi-layer fixed-thickness inversion.
### Release-20160428
	- Initial public release.

## Currently included programs
1. GAFORWARDMODELTDEM - 1D forward modelling program for AEM data
2. GALEISBSTDEM - deterministic 1D sample by sample inversion of AEM data
3. GARJMCMCTDEM - stochastic 1D sample by sample inversion of AEM data
4. CTLINEDATA2SGRID - (undocumented) convert inversion outputs to GoCAD SGrids 
5. CTLINEDATA2GEOREFIMAGE - (undocumented) convert inversion outputs to static georeferenced section images that can be displayed in a 2D GIS (a poor man's 3D).
6. CTLINEDATA2SLICEGRIDS - (undocumented) convert inversion outputs to layer, depth and elevation-slice grids in ErMapper format
7. CTLINEDATA2CURTAINIMAGE - (undocumented) convert inversion outputs to GA's Earth Sci curtain image format
8. EXAMPLE_FORWARD_MODEL - (undocumented) simple C++ example of how to use the code to run a forward models.


## Documentation
- [User Manual](docs/GA-AEM_Programs_User_Manual.pdf)
- [Theoretical details for GALEISBSTDEM](docs/GALEISBSTDEM_Inversion_Algorithm_Theoretical_Details.pdf)

## Building on linux

- Note that Windows 10 users can build and run the code easily on the *[Ubuntu 20.04](https://www.microsoft.com/en-au/p/ubuntu-2004/9n6svws3rx71#activetab=pivot:overviewtab)* emulator app.

### 1. Prerequisites

- You are likely to need to install the following or later packages on your linux system using "sudo apt install"

```bash
> sudo apt install cmake
> sudo apt install gcc-9
> sudo apt install libfftw3-dev
> sudo apt install libnetcdf-c++4-dev
> sudo apt install libopenmpi-dev
> sudo apt install libgdal-dev
> sudo apt install libpetsc-real-dev
```

### 2. Building using cmake  

1. The cmake program uses the file [*CMakeLists.txt*](CMakeLists.txt) to build the executables and libraries.
1. The simplest form of cmake
```bash
> cd [directory where ga-aem repo is]  
> mkdir my-build-dir                            // my-build-dir is a temporary directory for building
> cd my-build-dir                               // change to the build-dir
> cmake ..                                      // configure and generate using ../CMakeLists.txt in the directory above
> cmake --build all                             // build all targets
> cmake --install . --prefix my-install-dir     // install the executables, libraries and headers in my-install-dir 
```
1. To choose a specific compiler
	- Replace  the line 
```bash
> cmake ..
```
	- with, for example, one of these lines
```bash
> cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
> cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
```

1. Other hints

- See [*cmake_build_script_examples.sh*](cmake_build_script_examples.sh) for example
- You may disable PETSc, MPI, GDAL, NETCDF by using the following cmake command line, if for example, you do not have them installed or want to use them. 
```bash
> cmake -DUSE_PETSC=off -DUSE_MPI=off -DUSE_GDAL=off -DUSE_NETCDF=off etc etc ..
```

### 3. Building using makefiles

- Edit the file [makefiles/run_make.sh](makefiles/run_make.sh) to setup for your particular system and compiler
	- set the C compiler (e.g., cc=gcc or cc=icc)
	- set the C++ compiler (e.g., cxx=g++ or cxx=icpc)
	- set the MPI C++ compiler (e.g., mpicxx=mpiCC)
	- set the C++ compiler flags (e.g. cxxflags='-std=c++17 -O3 -Wall')
	- set the executable directory (e.g., exedir='../bin/gnu')
- Then to build 
```bash
> cd makefiles                              // change to makefile directory 
> run_make.sh gnu allclean                  // build with gnu compiler with allclean 
```
- Matlab shared library should go into ga-aem/matlab/bin/linux (.dll on Windows or .so on linux)
- Python shared library should go into ga-aem/python/gatdaem1d (.dll on Windows or .so on linux)

## Building on Windows
- You can build the programs with the free Microsoft Visual Studio Community 2019.
- Visual Studio solution and project files are supplied.
- Open ga-aem\visualstudio\ga-aem-all\ga-aem-all.sln to compile all programs plus the matlab and python shared libraries.
- Alternatively open individual program solutions files in their respective directories.

## Additional source code and library dependencies

The programs required additional code or libraries as described below.

1. CPP-UTILS
	- CPP-UTILS is a repository of C++ utility classes and functions that are used across this and several other other projects.
	- Up until Release-20160606 these source files were included in the src\ directory but were moved into a separate submodule so that they can be used and maintained easily across several different projects.
	- CPP-UTILS is included as a git submodule of this repository (see [submodules](submodules/README.md)).
	- Only required if you are compiling the code.
	- Not required if you are just going to use the precompiled executables.

2. Eigen 
	- A C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
	- Language: C++
	- Website: http://eigen.tuxfamily.org/.	
	- Repository: https://github.com/eigenteam/eigen-git-mirror.git	
	- Eigen is not required if you just want to run the precompiled Windows executables.
	- Eigen is only required if you want to compile the source code.
	- Eigen is included as a git submodule of this repository (see [submodules](submodules/README.md)).


3. Message Passing Interface (MPI)
	- MPI is a standard used for parallel computation.  The MPI standard has been implemented in several different flavours by different consortia.
	- See https://www.mpi-forum.org
	- Both Windows and Linux users will need to install a suitable MPI for your system.
	- MPI is used in the following parallel enabled programs:
		- galeisbstdem.exe
		- galeisbsfdem.exe
		- garjmcmctdem.exe
		- galeiallatonce.exe
	- The pre-compiled Windows inversion programs require Microsoft HPC Pack 2012 because that is the flavour of MPI they have been compiled and linked with.
	- They can be recompiled using MPICH or other flavours of MPI if required.
	- If you do not want to install MPI on your Windows system, you can use the galeisbstdem-nompi.exe, which has not been linked with MPI but uses OpenMP shared memory parallelism instead (see manual for details).

4. The Fastest Fourier Transform in the West (FFTW)
	- FFTW is an optimised Fast Fourier Transform package developed at MIT by Matteo Frigo and Steven Johnson.
	- See http://www.fftw.org
	- FFTW is required if you are compiling the time-domain forward modelling or inversion programs on Windows or Linux.
	- Linux users will need to install a suitable FFTW for your system.
	- Windows users can obtain the header files, precompiled libraries and the runtime dlls from GitHub here,
		- https://github.com/rcb547/fftw3.2.2.dlls.git, or
		- git clone git@github.com:rcb547/fftw3.2.2.dlls.git.
	- FFTW is required if you want to execute the precompiled time-domain forward modelling or inversion programs on Windows.
	- FFTW is also required if you are just going to use the precompiled executables on Windows.
	- The directory containing the 64-bit FFTW dlls need to be in your Windows search path.

5. Portable, Extensible Toolkit for Scientific Computation (PETSc)
	- PETSc, pronounced PET-see (the S is silent), is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.
	- See https://www.mcs.anl.gov/petsc.
	- PETSc is only used by the program galeiallatonce.exe.
	- PETSc is only required for compiling galeiallatonce.exe, not executing it.
	- PETSc is not required for executing the precompiled Windows programs.
	- Linux users will need to install PETSc on your system.
	- Windows users can obtain the header files and precompiled libraries from GitHub, here
		- https://github.com/rcb547/petsc-3.9.4-vs2017-win64-libraries, or
		- git clone git@github.com:rcb547/petsc-3.9.4-vs2017-win64-libraries.git.

