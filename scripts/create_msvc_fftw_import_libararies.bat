@ECHO OFF

REM Run this from the Visual Studio command prompt to get the correct paths before compiling with CMake on Visual Studio
SET FFTW_DIR=%LocalAppData%\fftw-3.3.5-dll64
CD %FFTW_DIR%
lib /machine:x64 /def:libfftw3-3.def
pause

