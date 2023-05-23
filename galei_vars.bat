@ECHO ---------------------------------------------
@ECHO Adding GA-LEI directories to PATH env
@ECHO ---------------------------------------------

@ECHO OFF
REM Add executable and library directories to your search path
set path=.\bin\x64\Release\;%path%
set path=C:\Program Files\fftw-3.2.2-dll64;%path%
set path=C:\Program Files\netCDF 4.9.2\bin;%path%
