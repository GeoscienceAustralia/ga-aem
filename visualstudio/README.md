# Microsoft Visual Studio projects for GA-AEM

## Description
On Windows systems you can build the programs with the [Microsoft Visual Studio](https://visualstudio.microsoft.com) 2019 (or later) software in a GUI based integrated development environment (IDE).  This is in fact how the code has been developed. 
- Be certain to select and install the ***`Desktop development with C++`*** workload in the Visual Studio installer.
- The visualstudio sub-directory contains the project, solution, and property-sheet files for building the various programs and libraries with the Microsoft Visual Studio C++ compiler (MSVC).
- However some path updates will be required depending on how/where you have installed the third-party dependency libraries.
- The programs are known to be successfully built with VS2019 and VS2022 on Windows 11.
- This directory has nothing to do with building using CMAKE.  Although the Visual Studio compiler can be used by CMAKE, the files in this directory are not used by CMAKE.  This directory is only for compiling in the Visual Studio IDE GUI program.
- These notes only apply to Windows systems.

## Contents
The sub-directories are organized as follows:
- bin\
	- Compiled and linked 64 bit program executables.
	- Release versions end up in x64\Release.
	- Debug versions end up in x64\Debug.
	- This is not the directory where the CMAKE builds end up.
- lib\
	- Compiled and linked 64 bit libraries.
	- Release versions end up in x64\Release.
	- Release versions end up in x64\Debug.
- propertysheets\
	- Property sheets that help set up of the compiler paths and options for various packages. Some path updates will be required depending on how/where you have installed the third-party dependency libraries on your system.
- ga-aem-all\
	- directory containing the MSVC solution file for building all projects (programs/libraries) with one solution. 

The following sub-folders contain the individual program/library project and solution files:
- ctlinedata2curtainimage\
- ctlinedata2georefimage\
- ctlinedata2sgrid\
- ctlinedata2slicegrids\
- example_forwardmodel_c\
- example_forward_model\
- gaforwardmodeltdem\
- galeiallatonce\
- galeisbsfdem\
- galeisbstdem\
- galeisbstdem-nompi\
- garjmcmctdem\
- gatdaem1d_c_library\
- gatdaem1d_matlab\
- gatdaem1d_python\
- libnetcdfcxx\
- removelog10conductivityfromsgrid\

## Build Instructions
To build all programs plus the matlab and python shared libraries.
1. Open Microsoft Visual Studio 2019 or later
2. File | Open Project
3. In the Open Project/Solution dialog box navigate to, select, and then open <ga-aem-repository-directory>\visualstudio\ga-aem-all\ga-aem-all.sln
4. Choose Build | Build Solution

Alternatively to build just an individual program, for example galeisbstdem
1. Open Microsoft Visual Studio 2019 or later
2. File | Open Project
3. In the Open Project/Solution dialog box navigate to, select, and then open <ga-aem-repository-directory>\visualstudio\galeisbstdem\galeisbstdem.sln
4. Choose Build | Build Solution

### Note on the Matlab and Python shared library outputs
- The Matlab shared library MEX file (.mexe64) will be output to the directory
	- <ga-aem-repository-directory>\matlab\bin\gatdaem1d.mexw64
- The Python shared library (.dll) file will be output to the directory
	- <ga-aem-repository-directory>\python\gatdaem1d\gatdaem1d.dll

### Environment variables
- The Visual Studio project files rely on some environment variable to be set so that the include and library files of the external packages can be found.  
- To build all programs the following variables would need to be set: FFTW_DIR, MATLAB_DIR, MSMPI_BIN, MSMPI_INC, MSMPI_LIB64, GDAL_DIR, GDAL_DATA, PROJ_LIB, NETCDF_DIR, PETSC_DIR, PETSC_BIN.  
- For example, on the Windows 11 syeditstem used for development, the following variables were set:
```text
	FFTW_DIR=%LocalAppData%\fftw-3.3.5-dll64
	MSMPI_BIN=C:\Program Files\Microsoft MPI\Bin\
	MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\
	MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\
	GDAL_DIR=%LocalAppData%\gdal-3.7.1-mapserver-8-0-1
	GDAL_DATA=%GDAL_DIR%\bin\gdal-data
	PROJ_LIB=%GDAL_DIR%\bin\proj9\SHARE
	NETCDF_DIR=C:\Program Files\netCDF 4.9.2
	PETSC_DIR=%LocalAppData%\petsc\3.9.4\vs2017
	PETSC_BIN=%PETSC_DIR%\win64_release\lib
	MATLAB_DIR=C:\Program Files\MATLAB\R2022b
```
- Typically %LocalAppData% is located at `C:\Users\<username>\AppData\Local`.
- These paths need to be adjusted according to where your particular packages are installed.

## Install Instructions
It might be convenient to run the script ***create_package.bat*** after building the executables and libraries.  
- This copies the various built executables, libraries and the examples folder to a standalone install directory.
- It is useful for distributing the directory to a work colleague for independent use.
- This is similar to, but independent of, the CMAKE install step.
- Change the variable INSTALL_DIR at the top of the batch file to a suitable directory path.
```bat
	SET INSTALL_DIR=%LocalAppData%\GA-AEM
	or
	SET INSTALL_DIR=C:\win10dev\GA-AEM
```