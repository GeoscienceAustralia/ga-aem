@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM If you do not have MPI installed substitute galeisbstdem-nompi.exe in place of galeisbstdem.exe

REM galeisbstdem-nompi.exe galeisbstdem-do-not-solve-geometry.con
REM galeisbstdem-nompi.exe galeisbstdem-solve-rxpitch.con
galeisbstdem-nompi.exe galeisbstdem-solve-rxpitch-and-offsets.con

REM Standalone on a single CPU
REM galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con
REM galeisbstdem.exe galeisbstdem-solve-rxpitch.con
REM galeisbstdem.exe galeisbstdem-solve-rxpitch-and-offsets.con

REM Use 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con 4
REM galeisbstdem.exe galeisbstdem-solve-rxpitch.con 4
REM galeisbstdem-nompi.exe galeisbstdem-solve-rxpitch-and-offsets.con 4

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-do-not-solve-geometry.con
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-solve-rxpitch.con
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-solve-rxpitch-and-offsets.con

pause


