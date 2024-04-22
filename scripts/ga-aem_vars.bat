@ECHO ===========================================================================
@ECHO Adding GA-AEM programs and dependencies to the PATH environment variable
@ECHO ===========================================================================

@ECHO OFF
REM Add executable and library directories to your search path
REM Modify the paths below to suit your computer's folder structure
REM Ideally you would set GA-AEM_ROOT as an environment variable in your user environment (e.g. Start | Edit environment variable for your account). After setting, be sure to open a fresh command window or Windows Explorer window.
set GA-AEM_ROOT=%LocalAppData%\GA-AEM

REM Note that GDAL will precede NETCDF in the PATH variable to prevent possible conflicts in the dlls being used by GDAL and NETCDF

REM The GA-AEM compiled programs and Matlab and Python interface dlls
set GA-AEM_PATH=%GA-AEM_ROOT%\bin;%GA-AEM_ROOT%\matlab\bin;%GA-AEM_ROOT%\python\gatdaem1d

REM FFTW NOT required for the ctlinedata* programs
set GA-AEM_PATH=%GA-AEM_PATH%;%LocalAppData%\fftw-3.3.5-dll64

REM GDAL only required for ctlinedata2slicegrids.exe and ctlinedata2curtainimage.exe
REM PROJ_LIB is in-turn used by GDAL for geographic datum/projection functionality
set GA-AEM_PATH=%GA-AEM_PATH%;%LocalAppData%\gdal-3.7.1-mapserver-8-0-1\bin
set PROJ_LIB=%LocalAppData%\gdal-3.7.1-mapserver-8-0-1\bin\proj9\share

REM NETCDF required by garjmcmctdem.exe, and only optionally required for galeisbstdem.exe, galeisbstdem-nompi.exe
set GA-AEM_PATH=%GA-AEM_PATH%;C:\Program Files\netCDF 4.9.2\bin

REM PETSc only required for galeiallatonce.exe
set GA-AEM_PATH=%GA-AEM_PATH%;%LocalAppData%\petsc\3.9.4\vs2017\win64_release\lib

REM Finally prepend %GA-AEM_PATH% to the existing %PATH% to prevent conflicts with other version installed by say Python or QGIS/ArcGIS
set PATH=%GA-AEM_PATH%;%PATH%

REM Additional while not PATH variables, GDAL_DATA and PROJ_LIB may also be required at runtime by GDAL for geographic datum/projection functionality
set GDAL_DATA=%LocalAppData%\gdal-3.7.1-mapserver-8-0-1\bin\gdal-data
set PROJ_LIB=%LocalAppData%\gdal-3.7.1-mapserver-8-0-1\bin\proj9\share

REM ECHO %GA-AEM_ROOT%
REM ECHO %GA-AEM_PATH%
REM ECHO %PATH% > path_variable.txt

REM These variables below should only be required for the CMake or Visual Studio build process rather than running programs. However the PATH variables need to be set (as per above and suitably modified for your system) on a system where the programs are being executed.
REM set FFTW_DIR=%LocalAppData%\fftw-3.3.5-dll64
REM set GDAL_DIR=%LocalAppData%\gdal-3.7.1-mapserver-8-0-1
REM set GDAL_DATA=%GDAL_DIR%\bin\gdal-data
REM set PROJ_LIB=%GDAL_DIR%\bin\proj9\SHARE
REM set NETCDF_DIR=C:\Program Files\netCDF 4.9.2
REM set PETSC_DIR=%LocalAppData%\petsc\3.9.4\vs2017
REM set PETSC_BIN=%PETSC_DIR%\win64_release\lib
REM set MATLAB_DIR=C:\Program Files\MATLAB\R2022
REM set MSMPI_BIN=C:\Program Files\Microsoft MPI\Bin
REM set MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\
REM set MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\

