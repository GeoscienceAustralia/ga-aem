@ECHO ---------------------------------------------------------------------------
@ECHO Adding GA-AEM programs and dependencies to the PATH environment variable
@ECHO ---------------------------------------------------------------------------

@ECHO OFF
REM Add executable and library directories to your search path
REM Modify the paths below to suit your computer's folder structure
REM Ideally you would set GA-AEM_ROOT and the PATH variable in your user environment (e.g. Start | Edit environment variable for your account)

set path=%GA-AEM_ROOT%\bin;%PATH%

REM FFTW NOT required for the ctlinedata* programs
set path=%LocalAppData%\fftw-3.3.5-dll64;%PATH%

REM NetCDF only required for galeisbstdem.exe, galeisbstdem-nompi.exe, garjmcmctdem.exe
set path=C:\Program Files\netCDF 4.9.2\bin;%PATH%

REM GDAL only required for ctlinedata2georefimage.exe and ctlinedata2curtainimage.exe
set path=%LocalAppData%\gdal-3.0.4\bin;%PATH%
REM set path=%LocalAppData%\%LocalAppData%\gdal-3.7.1-mapserver-8-0-1;%PATH%

REM PETSc only required for galeiallatonce.exe
SET path=%LocalAppData%\petsc\3.9.4\vs2017\win64_release\lib;%PATH%

