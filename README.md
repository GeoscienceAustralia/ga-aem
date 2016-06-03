# GA-AEM Source Code Repository

## Geoscience Australia Airborne Electromagnetics Programs

- Author:	Ross C Brodie, Geoscience Australia (ross.c.brodie at ga.gov.au)
- Language:	mostly C++, some matlab, some python

## Releases
### Release-20160603
	- Added Python 3 interface for simple forward modelling and derivatives only.
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
- User manual (see ga-aem/docs/GA AEM Programs User Manual.pdf)
- Theoretical details for GALEISBSTDEM (see ga-aem/docs/GALEISBSTDEM Inversion Algorithm Theoretical Details .pdf)

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

## Third party software dependencies
1. Template Numerical Toolkit (TNT)
	- see http://math.nist.gov/tnt/index.html
	- TNT is included in this repositiry, there is no need to download it.

2. The Fastest Fourier Transform in the West (FFTW)
	- see http://www.fftw.org
	- Windows binaries (dlls) are include in the repository.
	- Linux users will need to install a suitable FFTW for your system.

3. Message Passing Interface (MPI)
	- see https://www.mpi-forum.org
	- Both Windows and Linux users will need to install a suitable MPI for your system.
	- The pre-compiled Windows inversion programs require Microsoft HPC Pack 2012 because thats the flavour of MPI they have been compiled and linked with.
		- They can be recompiled using MPICH or other flouvours of MPI if required.
		- If you do not want to install MPI on your Windows system, you can use the galeisbstdem-nompi.exe.
		- This program has not been linked with MPI but uses OpenMP shared memory parallelism instead (see manual for details).

