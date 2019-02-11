@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)

REM set path=C:\Microsoft_HPC_Pack_2012\Bin;%path%
REM set path=..\..\..\..\fftw3.2.2.dlls\64bit;%path%
set path=..\..\..\..\petsc\3.9.4\vs2017\win64_release\lib;%path%
set path=..\..\..\bin\x64\Release\;%path%

REM Standalone
REM #galeisbstdem.exe galeisbstdem.con

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem.con

REM Use 4 OpenMP threads
galeisbstdem.exe galeisbstdem.con 4

REM If you don't have MPI installed use the non-MPI version
REM Use 4 OpenMP threads 
REM galeisbstdem-nompi.exe galeisbstdem.con 4

pause


