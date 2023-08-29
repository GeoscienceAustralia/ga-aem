@ECHO ---------------------------------------------------------------------------
@ECHO Adding GA-AEM programs and dependencies to the PATH environment variable
@ECHO ---------------------------------------------------------------------------

@ECHO OFF
REM Add executable and library directories to your search path
REM Modify the paths below to suit your computer's folder structure
set path=%GA-AEM_ROOT%\bin;%PATH%
set path=C:\Program Files\netCDF 4.6.1\bin;%PATH%
set path=%LocalAppData%\fftw-3.3.5-dll64;%path%
set path=%LocalAppData%\gdal-3.0.4\bin;%PATH%
rem set path=%LocalAppData%\gdal-3.0.4\bin\gdal\plugins;%path%
