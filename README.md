# GA-AEM Source Code Repository

# Description
This is a repository for Geoscience Australia's programs and utilities for forward modelling and inversion of Airborne Electromagnetic (AEM) data. It includes Matlab and Python interfaces for forward model and derivative calculations. It also includes some programs for post-processing inversion results into sections, images and grids.

## Authors
- Ross C Brodie, Geoscience Australia (formerly)
- Richard Taylor, Geoscience Australia (formerly)

## Languages
- Mostly C++.
- Some Matlab.
- Some Python.

## Included content

### User programs
- gaforwardmodeltdem.exe - 1D forward modelling program for time-domain AEM data
- galeisbstdem.exe - deterministic 1D sample-by-sample or bunch-by-bunch inversion of time-domain AEM data
- galeisbstdem-nompi.exe - as above but without any parallel MPI support or dependency
- garjmcmctdem.exe - (undocumented) stochastic 1D sample by sample inversion of time-domain AEM data
- galeiallatonce.exe - (undocumented) stochastic 1D sample by sample inversion of time-domain AEM data
- ctlinedata2sgrid.exe - (undocumented) convert inversion outputs to GoCAD SGrids 
- ctlinedata2georefimage.exe - (undocumented) convert inversion outputs to static georeferenced section images that can be displayed in a 2D GIS (a poor man's 3D).
- ctlinedata2slicegrids.exe - (undocumented) convert inversion outputs to layer, depth and elevation-slice grids in ErMapper format
- ctlinedata2curtainimage.exe - (undocumented) convert inversion outputs to GA's Earth Sci curtain image format
- removelog10conductivityfromsgrid.exe - legacy program for removing the log10 conductivity from GoCAD SGrids
### User examples
- Examples of how to use the programs for various AEM systems
### For Matlab users
- Matlab interface via MEX file (shared library) with examples
- See [*here*](matlab/README.md) for details.
### For Python users
- Python interface via shared library with examples.
- See [*here*](python/README.md) for details.
### For developers/coders
- example_forward_model.exe - (for developers) simple C++ language example of how to use the code in C++ to run a forward models.
- example_forward_model_c.exe - (for developers) simple C language example of how to use the code in C++ to run a forward models.

### User documentation
- [User Manual](docs/GA-AEM_Programs_User_Manual.pdf)
- [Theoretical details for GALEISBSTDEM](docs/GALEISBSTDEM_Inversion_Algorithm_Theoretical_Details.pdf)

# Cloning the repository
When initially cloning the repository in git you should use the `--recursive` option so that all of the submodules and their respective submodules are initialised and populated with code.
```bash
> git clone --recursive https://github.com/GeoscienceAustralia/ga-aem.git
```
or if you use SSH authentification,
```bash
> git clone --recursive git@github.com:GeoscienceAustralia/ga-aem.git
```

## Submodules
The ga-aem project has several source code dependencies that are included as git submodules from other open-source projects. The submodules are only required for building the programs and are not required if you are just using precompiled executables. See [*here*](./submodules/README.md) for details of how the submodules should be initialised and updated.

## Third-party library dependencies
For full functionality and to build all programs the following packages are required: FFTW, MPI, NetCDF, GDAL and PETSc.
- See [*here*](README-Dependencies.md) for details of how to obtain and install the library dependencies.
- Not all the dependencies are required for all the programs, as detailed below,
	- FFTW
		- required for galeisbstdem.exe, galeisbstdem-nompi.exe, garjmcmctdem.exe, galeiallatonce.exe, Matlab and Python interfaces.
	- MPI
		- optional for galeisbstdem.exe
		- required for garjmcmctdem.exe and galeiallatonce.exe
	- NetCDF
		- optional for galeisbstdem.exe and galeisbstdem-nompi.exe
		- required for garjmcmctdem.exe
	- GDAL
		- required only for ctlinedata2slicegrids.exe and ctlinedata2curtainimage.exe
	- PETSc
		- required only for galeiallatonce.exe

# Building (compiling and linking) the code
- The programs can be built from source code on both Linux and Windows systems, and probably other architectures.
- The build system uses CMake (>=v3.16) software.
- Traditional Makefiles are deprecated in this release of ga-aem.
- It is typically simpler to build the code on Linux, however note that Windows users can build and run the code easily on the free [*Ubuntu 20.04*](https://www.microsoft.com/en-au/p/ubuntu-2004/9n6svws3rx71#activetab=pivot:overviewtab) emulator app available free from the Microsoft Store.
- Nevertheless, the code definitely can be built on Windows with CMake or with the Microsoft Visual Studio IDE.

## Building using CMake
- CMake can be downloaded from https://cmake.org/download.
- The CMake program uses the file [*CMakeLists.txt*](CMakeLists.txt) to build the executables and libraries.  Unless you really know what you are doing do not edit this file.
- CMake involves a generate step, a build step, and an install step.
- The most basic way to use CMake is as follows:
```bash
> cd [directory where ga-aem repository is located]
> mkdir <build-dir>                                // <build-dir> is a temporary directory for building
> cd <build-dir>                                   // change to the <build-dir>
> cmake -DCMAKE_BUILD_TYPE=Release ..              // generate the cache using ../CMakeLists.txt in the directory above
> cmake --build .                                  // build all targets
> cmake --install . --prefix <install-dir>         // install the executables, libraries, headers, docs and Matlab and Python into <install-dir>
```
- The above should work on Linux if you have a C and C++ compiler installed along with all the dependencies.
- However, it is recommend that you inspect, copy and then alter one of the four provided CMake build script examples below to suit your specific purposes.
	- MSVC on Windows [*cmake_build_script_windows-vs2022.bat*](cmake_build_script_windows-vs2022.bat)
	- GNU compiler on Ubuntu (including emulator on Windows) [*cmake_build_script_ubuntu.sh*](cmake_build_script_ubuntu.sh)
	- GNU compiler on Gadi cluster [*cmake_build_script_gadi-gnu.sh*](cmake_build_script_gadi-gnu.sh)
	- Intel compiler on Gadi cluster [*cmake_build_script_gadi-intel.sh*](cmake_build_script_gadi-intel.sh)
### CMake generate step
- To choose a specific compiler, replace the line,
	```bash
	> cmake -DCMAKE_BUILD_TYPE=Release ..
	```
	with, for example, one of these lines for the GNU, Intel and Intel-LLVM compilers respectively,
	```bash
	> cmake -DCMake_C_COMPILER=gcc -DCMake_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ..
	> cmake -DCMake_C_COMPILER=icc -DCMake_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release ..
	> cmake -DCMake_C_COMPILER=icx -DCMake_CXX_COMPILER=icpx -DCMAKE_BUILD_TYPE=Release ..
	```
	or for Microsoft Visual Studio 2019 or 2022 compilers respectively,
	```dos
	> cmake -G "Visual Studio 16 2019" -A x64 -DCMAKE_CXX_COMPILER=msvc -DCMAKE_BUILD_TYPE=Release ..
	> cmake -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_COMPILER=msvc -DCMAKE_BUILD_TYPE=Release ..
	```
- Specific dependencies may be disabled by using one or more of the following optional command line switches in the initial generate step. These options are ON by default.
```text
	-DWITH_FFTW=OFF
	-DWITH_MPI=OFF
	-DWITH_GDAL=OFF
	-DWITH_NETCDF=OFF
	-DWITH_PETSC=OFF
```
This may be useful if, for example, you do not have them installed or do not need that particular functionality. For example a minimalistic version of `galeisbstdem.exe` could be built with the following in the generate step,
```bash
	cmake -G "Visual Studio 17 2022" -A x64 -Wno-dev -DCMAKE_CXX_COMPILER=msvc -DWITH_FFTW=ON -DWITH_MPI=OFF -DWITH_GDAL=OFF -DWITH_NETCDF=OFF -DWITH_PETSC=OFF -DCMAKE_BUILD_TYPE=Release ..
```
### CMake build step
- In the build step all targets (executables/libraries) can be built,
	```bash
	> cmake --build . --config=Release
	```
	or only specific targets can be specified,
	```bash
	> cmake --build . --target galeisbstdem --config=Release
	> cmake --build . --target ctlinedata2sgrid --config=Release
	```
- Note that in the build step, the *--config=Release* switch is required for the Visual Studio (and other [*generator*](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html) style) compilers, but it is not required for the GNU and Intel compilers on Linux.
### CMake Install step
- The install step installs (copies) all the executables, libraries, headers, examples, Matlab and Python interfaces to an installation directory.
	```bash
	cmake --install . --prefix <install-dir>
	```
- The <install-dir> is a suitable installation directory, for example on Windows it might be,
	```bash
	cmake --install . --prefix %LocalAppData%\GA-AEM
	```
## Building on Windows with the Microsoft Visual Studio IDE
- On Windows systems you can build the programs with the free Microsoft Visual Studio 2019 (or later) software in a GUI based integrated development environment (IDE).  This is in fact how the code has been developed.
- For convenience Microsoft Visual Studio project, solution and property sheet files are supplied. However some path updates will be required depending on how/where you have installed the third-party dependency libraries.
- See [*here*](visualstudio/README.md) for more details.

# Releases

## Release-20240424
- To be completed

## Release-20160606
- Added Python 3.x interface for simple forward modelling and derivatives only.
- Added Matlab interface for simple forward modelling and derivatives only.
- Changed how the PPM normalisation is carried out. Now PPM normalisation is by directional-component-wise with respect to the maximum primary dB/dt or B-field at the receiver for a reference system geometry (which is usually estimated on a per flight or per survey basis). Previously PPM normalisation was with respect to the system geometry for the forward model being run.
- Added GEOTEM (1996 ppm system) and SPECTREM (ppm system) examples.
- Fixed a bug in the thickness derivative of the second bottom layer. This may have effected few-layer inversions, but not multi-layer fixed-thickness inversion.
## Release-20160428
- Initial public release.


