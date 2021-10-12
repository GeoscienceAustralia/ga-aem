@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=C:\Microsoft_HPC_Pack_2012\Bin;%path%
REM set path=..\..\..\..\fftw3.2.2.dlls\64bit;%path%
set path=..\..\..\..\petsc\3.9.4\vs2017\win64_release\lib;%path%
set path=..\..\..\bin\x64\Release\;%path%

REM Using 4 MPI processes
mpiexec -np 4 galeiallatonce.exe galeiallatonce.con

pause

