@ECHO off

REM Run this script from the "X64 Native Tools Command Prompt for VS 2019 (or VS 2022)" so that all the compiler stuff is set up

rem BUILD_DIR is a temporary directory for building (compiling and linking)
set BUILD_DIR=%cd%\build-windows-vs2022
rem INSTALL_DIR is the directory for installing the built package and examples (e.g. c:\myprograms\ga-aem)
set INSTALL_DIR=%cd%\install-windows-vs2022

rem Set the FFTW path if necesary
SET FFTW_DIR=C:\Users\rossc\work\code\third_party\fftw3.2.2.dlls\64bit

rem Set the PETSc path if necesary
SET PETSC_DIR=C:\Users\rossc\work\code\third_party\petsc-installed\3.9.4\vs2017
SET PETSC_INCLUDE_DIRS=%PETSC_DIR%\include;%PETSC_DIR%\win64_release\include
SET PETSC_LIB_DIR=%PETSC_DIR%\win64_release\lib
SET PETSC_LIBRARIES=libf2cblas.lib;libf2clapack.lib;libpetsc.lib
SET PETSC_LDFLAGS=

rem Set the NetCDF path if necesary
SET NETCDF_DIR=C:\Program Files\netCDF 4.6.1
SET NETCDF_INCLUDE_DIR=%NETCDF_DIR%\include
SET NETCDF_LIB_DIR=%NETCDF_DIR%\lib
SET NETCDF_LIBRARIES=netcdf.lib;hdf5.lib;hdf5_hl.lib;libcurl_imp.lib

rem Maybe delete the BUILD_DIR for a fresh start
rem RMDIR /S /Q %BUILD_DIR%

rem Create and cd to the BUILD_DIR
mkdir %BUILD_DIR%
cd    %BUILD_DIR%

rem First generate the build cache first
rem cmake -Wno-dev -G "Visual Studio 16 2019" -A x64 -DCMAKE_CXX_COMPILER=msvc ..
cmake -Wno-dev -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_COMPILER=msvc ..

rem Switches for turning off certain dependencies if they are not wanted or available
rem	-DWITH_MPI=OFF
rem	-DWITH_NETCDF=OFF
rem	-DWITH_GDAL=OFF
rem	-DWITH_PETSC=OFF
rem cmake -Wno-dev -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_COMPILER=msvc -DWITH_MPI=OFF -DWITH_NETCDF=OFF -DWITH_GDAL=OFF -DWITH_PETSC=OFF ..

rem Build and install everything
cmake --build . --config=Release
cmake --install . --prefix %INSTALL_DIR%

rem Or alternatively ...

rem Build only particular targets
rem cmake --build . --target galeisbstdem --config=Release
rem cmake --build . --target garjmcmctdem --config=Release
rem cmake --build . --target galeisbsfdem --config=Release
rem cmake --build . --target galeiallatonce --config=Release

rem cmake --build . --target gaforwardmodeltdem --config=Release
rem cmake --build . --target example_forward_model --config=Release
rem cmake --build . --target example_forward_model_c --config=Release
rem cmake --build . --target gatdaem1d-static --config=Release

rem gatdaem1d-shared is required FOR python and matlab interfaces
rem cmake --build . --target gatdaem1d-shared --config=Release  

rem cmake --build . --target ctlinedata2sgrid --config=Release
rem cmake --build . --target ctlinedata2slicegrids --config=Release
rem cmake --build . --target ctlinedata2georefimage --config=Release
rem cmake --build . --target ctlinedata2curtainimage --config=Release
rem cmake --build . --target removelog10conductivityfromsgrid --config=Release

rem cmake --install . --prefix %INSTALL_DIR%

cd ..
PAUSE

