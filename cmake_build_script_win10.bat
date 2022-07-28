@ECHO off

rem %CommSpec% /k "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
SET Boost_ROOT=C:\Users\rossc\work\code\third_party\boost_1_79_0
SET CGAL_ROOT=C:\Users\rossc\work\code\third_party\CGAL-5.4.1
SET FFTW_ROOT=C:\Users\rossc\work\code\third_party\fftw3.2.2.dlls\64bit
rem SET NETCDFCXX_ROOT=C:\Users\rossc\work\code\repos\ga-aem-csiro\submodules\netcdf-cxx4

rem set PATH=%PATH%;C:\ProgramData\msys64\mingw64\bin;C:\ProgramData\msys64\usr\bin;D:\Software\ga-aem\third_party\fftw3.2.2.dlls\64bit
rem set VCTargetsPath=C:\Program Files (x86)\MSBuild\Microsoft.Cpp\v4.0\V140

rem set PATH=%PATH%;C:\ProgramData\msys64\mingw64\bin;C:\ProgramData\msys64\usr\bin;D:\Software\ga-aem\third_party\fftw3.2.2.dlls\64bit
rem set VCTargetsPath=C:\Program Files (x86)\MSBuild\Microsoft.Cpp\v4.0\V140

set BUILD_DIR=%cd%\win10-build-vs2019
set INSTALL_DIR=%cd%\win10-install-vs2019

rem RMDIR /S /Q %BUILD_DIR%
rem mkdir %BUILD_DIR%

cd %BUILD_DIR%
cmake -Wno-dev -G "Visual Studio 16 2019" -A x64 -DCMAKE_CXX_COMPILER=msvc ..
rem cmake --build . --target libnetcdfcxx --config=Release
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

PAUSE

