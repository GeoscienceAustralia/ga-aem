@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%


REM Standalone
galeisbstdem.exe galeisbstdem.con

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem.con

REM Use 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem.con 4

REM If you don't have MPI installed use the non-MPI version
REM Use 4 OpenMP threads 
REM galeisbstdem-nompi.exe galeisbstdem.con 4

pause


