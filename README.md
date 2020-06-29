# GA-AEM Source Code Repository
[![Build Status](https://travis-ci.com/GeoscienceAustralia/ga-aem.svg?branch=develop)](https://travis-ci.com/GeoscienceAustralia/ga-aem)

## Geoscience Australia Airborne Electromagnetics Programs

- Author:	Ross C Brodie, Geoscience Australia (ross.c.brodie at ga.gov.au)
- Language:	mostly C++, some matlab, some python

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

## Documentation
- [User Manual](docs/GA-AEM_Programs_User_Manual.pdf)
- [Theoretical details for GALEISBSTDEM](docs/GALEISBSTDEM_Inversion_Algorithm_Theoretical_Details.pdf)

## Building on Linux
- cd makefiles
- edit the file run_make.sh to setup for your compiler
	- set the C++ compiler (e.g., cxx=g++)
	- set the MPI C++ compiler (e.g., mpicxx=mpiCC)
	- set the C++ compiler flags (e.g. cxxflags='-std=c++11 -O3 -Wall')
	- set the executable directory (e.g., exedir='../bin/raijin/gnu')
- run_make.sh
- Matlab shared library should go into ga-aem/matlab/bin/linux (.dll on Windows or .so on linux)
- Python shared library should go into ga-aem/python/gatdaem1d (.dll on Windows or .so on linux)

## Building on Windows
- You can build the programs with the free Microsoft Visual Studion 2013 Express.
- Visual Studio solution and project files are supplied.
- Open ga-aem\vs2013\ga-aem-all\ga-aem-all.sln to compile all programs plus the matlab and python shared libraries.
- Alternatively open individual program solutions files in their respective directories.

## Additional source code and library dependencies

The programs required additional code or libraries as described below.

1. CPP-UTILS
	- CPP-UTILS is a repository of C++ utility classes and functions that are used across this and several other other projects.
	- Up until Release-20160606 these source files were included in the src\ directory but were moved into a separate submodule so that they can be used and maintained easily across several different projects.
	- CPP-UTILS is included as a git submodule of this repository (see [submodules](submodules/README.md)).
	- Only required if you are compiling the code.
	- Not required if you are just going to use the precompiled executables.

2. Template Numerical Toolkit (TNT)
	- TNT is a C++ linear algebra package developed by the National Institute of Standards and Technology (NIST).
	- See http://math.nist.gov/tnt/index.html.
	- TNT is included as a git submodule of this repository (see [submodules](submodules/README.md)).
	- TNT is only required if you are compiling the code.
	- TNT is not required if you are just going to use the precompiled executables.

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
	- Windows users can obtain the header files and precompiled libraries from the GitHub, here
		- https://github.com/rcb547/petsc-3.4.3-vs2013, or
		- git clone git@github.com:rcb547/petsc-3.4.3-vs2013.git.

