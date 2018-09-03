@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM Use 2 MPI processes
mpiexec -np 2 galeiallatonce.exe galeiallatonce.con

pause


