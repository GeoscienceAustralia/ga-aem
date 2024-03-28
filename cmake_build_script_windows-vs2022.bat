@ECHO off

REM Run this script from the "X64 Native Tools Command Prompt for VS 2019 (or VS 2022)" so that all the compiler stuff is set up

REM INSTALL_DIR is the directory for installing the built package and examples (e.g. c:\myprograms\ga-aem)
set INSTALL_DIR=%LocalAppData%\GA-AEM

REM BUILD_DIR is a temporary directory for building (compiling and linking)
set BUILD_DIR=%cd%\build-windows-vs2022

REM Optionally delete the BUILD_DIR to ensure a clean cache/start
REM DEL /S /Q %BUILD_DIR%

REM Set the FFTW path if not set in user environment externally
SET FFTW_DIR=%LocalAppData%\fftw-3.3.5-dll64

REM Set the NetCDF path if not set in user environment externally
SET NETCDF_DIR=C:\Program Files\netCDF 4.9.2

REM Set the GDAL path if not set in user environment externally
SET GDAL_DIR=%LocalAppData%\gdal-3.0.4

REM Set the PETSc paths if not set in user environment externally
SET PETSC_DIR=%LocalAppData%\petsc\3.9.4\vs2017
SET PETSC_INCLUDE_DIRS=%PETSC_DIR%\include;%PETSC_DIR%\win64_release\include
SET PETSC_LIB_DIR=%PETSC_DIR%\win64_release\lib
SET PETSC_LIBRARIES=libf2cblas.lib;libf2clapack.lib;libpetsc.lib
SET PETSC_LDFLAGS=

REM Create and cd to the BUILD_DIR
mkdir %BUILD_DIR%
cd    %BUILD_DIR%

REM First generate the build cache first
REM cmake -Wno-dev -G "Visual Studio 16 2019" -A x64 -DCMAKE_CXX_COMPILER=msvc ..
cmake -Wno-dev -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_COMPILER=msvc -DCMAKE_BUILD_TYPE=Release ..

REM Switches for turning off certain dependencies if they are not wanted or available
REM 	-DWITH_MPI=OFF
REM 	-DWITH_NETCDF=OFF
REM 	-DWITH_GDAL=OFF
REM 	-DWITH_PETSC=OFF
REM cmake -Wno-dev -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_COMPILER=msvc -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF -DWITH_NETCDF=OFF -DWITH_GDAL=OFF -DWITH_PETSC=OFF ..

REM Build and install everything
cmake --build . --config=Release
cmake --install . --prefix %INSTALL_DIR%

REM Or alternatively ...

REM Build only particular targets
REM cmake --build . --target galeisbstdem --config=Release
REM cmake --build . --target galeisbstdem-nompi --config=Release
REM cmake --build . --target garjmcmctdem --config=Release
REM cmake --build . --target galeisbsfdem --config=Release
REM cmake --build . --target galeiallatonce --config=Release

REM cmake --build . --target gaforwardmodeltdem --config=Release
REM cmake --build . --target example_forward_model --config=Release
REM cmake --build . --target example_forward_model_c --config=Release
REM cmake --build . --target gatdaem1d-static --config=Release

REM gatdaem1d-shared is required FOR python and matlab interfaces
REM cmake --build . --target gatdaem1d-shared --config=Release

REM cmake --build . --target ctlinedata2sgrid --config=Release
REM cmake --build . --target ctlinedata2slicegrids --config=Release
REM cmake --build . --target ctlinedata2georefimage --config=Release
REM cmake --build . --target ctlinedata2curtainimage --config=Release
REM cmake --build . --target removelog10conductivityfromsgrid --config=Release

REM cmake --install . --prefix %INSTALL_DIR%

cd ..
PAUSE

