@ECHO off

rem %CommSpec% /k "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
rem SET Boost_DIR=C:\Users\rossc\work\code\third_party\boost_1_79_0
rem SET CGAL_DIR=C:\Users\rossc\work\code\third_party\CGAL-5.4.1
SET FFTW_DIR=C:\Users\rossc\work\code\third_party\fftw3.2.2.dlls\64bit
SET PETSC_DIR=C:\Users\rossc\work\code\third_party\petsc-installed\3.9.4\vs2017
SET PETSC_INCLUDE_DIRS=%PETSC_DIR%\include;%PETSC_DIR%\win64_release\include
SET PETSC_LIB_DIR=%PETSC_DIR%\win64_release\lib
SET PETSC_LIBRARIES=libf2cblas.lib;libf2clapack.lib;libpetsc.lib
SET PETSC_LDFLAGS=

SET NETCDF_DIR=C:\Program Files\netCDF 4.6.1
SET NETCDF_INCLUDE_DIR=%NETCDF_DIR%\include
SET NETCDF_LIB_DIR=%NETCDF_DIR%\lib
SET NETCDF_LIBRARIES=netcdf.lib;hdf5.lib;hdf5_hl.lib;libcurl_imp.lib

set BUILD_DIR=%cd%\win10-build-vs2019
set INSTALL_DIR=%cd%\win10-install-vs2019

rem RMDIR /S /Q %BUILD_DIR%
mkdir %BUILD_DIR%

cd %BUILD_DIR%
cmake -Wno-dev -G "Visual Studio 16 2019" -A x64 -DCMAKE_CXX_COMPILER=msvc ..
rem cmake --build . --target galeisbstdem --config=Release
rem cmake --build . --target garjmcmctdem --config=Release
rem cmake --build . --target galeisbsfdem --config=Release
cmake --build . --target galeiallatonce --config=Release

rem cmake --build . --target gaforwardmodeltdem --config=Release
rem cmake --build . --target example_forward_model --config=Release
rem cmake --build . --target example_forward_model_c --config=Release
rem cmake --build . --target gatdaem1d-static --config=Release
rem cmake --build . --target gatdaem1d-shared --config=Release

rem cmake --build . --target ctlinedata2georefimage --config=Release
rem cmake --build . --target ctlinedata2curtainimage --config=Release
rem cmake --build . --target ctlinedata2sgrid --config=Release
rem cmake --build . --target ctlinedata2slicegrids --config=Release
rem cmake --build . --target removelog10conductivityfromsgrid --config=Release

rem cmake --build . --config=Release
rem cmake --install . --prefix %INSTALL_DIR%

cd ..
PAUSE

