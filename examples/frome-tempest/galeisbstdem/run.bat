@echo off

REM Add executable and FFTW directories to your search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM If you do not have MPI installed substitute galeisbstdem-nompi.exe in place of galeisbstdem.exe

REM galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con
galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con
REM galeisbstdem.exe galeisbstdem-xzcomponents-do-not-solve-geometry.con

REM Standalone on a single CPU

REM Use 4 OpenMP threads

REM Use 4 MPI processes

pause


