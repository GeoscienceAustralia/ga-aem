# GA-AEM Source Code Repository
[![Build Status](https://travis-ci.com/GeoscienceAustralia/ga-aem.svg?branch=develop)](https://travis-ci.com/GeoscienceAustralia/ga-aem)

## Geoscience Australia Airborne Electromagnetics Programs

- Author:	Ross C Brodie, Geoscience Australia (ross.c.brodie at ga.gov.au) and Richard Taylor, Geoscience Australia
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
Building GA-AEM on Linux requires that you install the following software:
- [CMake](https://cmake.org/) >= 3.12
- [netCDF-4](https://www.unidata.ucar.edu/software/netcdf/) with C++ bindings and headers
- [fftw3](http://www.fftw.org/) with headers

For production-scale inversions, you will also need a working MPI installation on your system, and the MPI headers must be on your compiler's search path.

There are also legacy shell scripts and makefiles in the makefiles subdirectory for use with make. These will require you to manually specify which compiler you would like to use (GCC or Intel compiler) and are not actively maintained, so the CMake build system is recommended.

To run an out-of-source build, avoiding polluting the source tree, run the following commands from the root of the repository:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
sudo make install
```
This will configure and compile the package and install the binaries into `/usr/local/bin`. If you'd like to install somewhere else or you aren't a `sudo`er you can specify a custom install prefix using `-DCMAKE_INSTALL_PREFIX`. The `-DCMAKE_BUILD_TYPE=Release` flag sets your compiler to optimise the compiled binaries for speed (e.g. `-O3` with gcc). If you would like to include debugging symbols in the binaries, you can set `-DCMAKE_BUILD_TYPE=Debug` instead.

Finally, when invoking `cmake`, you can configure the build manually to not use MPI by setting `-DWITH_MPI=OFF`, and to use OpenMP multi-threading with the GARJMCMC parallel tempering program by setting `-DWITH_OMP_GARJMCMC=ON`. These options can be useful if you're building GA_AEM to do synthetic inversions or single soundings on a desktop machine instead of production inversions on a cluster.

The CMake build system should also be capable of generating cross-platform recipes to compile GA-AEM on Windows or MacOS, but this is not tested.

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
