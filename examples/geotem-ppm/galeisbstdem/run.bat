@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM If you don't have MPI installed use the non-MPI version
REM Use 4 OpenMP threads 
galeisbstdem-nompi.exe galeisbstdem_z_nosolvegeometry.con 4
rem galeisbstdem-nompi.exe galeisbstdem_xz_nosolvegeometry.con 4
//Solving for geometry on secondary field data is very "iffy" because it is the primary field that is more sensitive to geometry.
rem galeisbstdem-nompi.exe galeisbstdem_xz_solvegeometry.con 4


REM Standalone
rem galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con
rem galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con
rem galeisbstdem.exe galeisbstdem_xz_solvegeometry.con

REM Use 4 MPI processes
rem mpiexec -np 4 galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con
rem mpiexec -np 4 galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con
rem mpiexec -np 4 galeisbstdem.exe galeisbstdem_xz_solvegeometry.con

REM Use 4 OpenMP threads
rem galeisbstdem.exe galeisbstdem_z_nosolvegeometry.con 4
rem galeisbstdem.exe galeisbstdem_xz_nosolvegeometry.con 4
rem galeisbstdem.exe galeisbstdem_xz_solvegeometry.con 4

pause

