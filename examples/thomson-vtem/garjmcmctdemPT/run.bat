echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM Run as standalone process
REM garjmcmctdem.exe garjmcmctdem.con

REM Run with 4 MPI processes
REM mpiexec -np 4 garjmcmctdem.exe garjmcmctdem.con

garjmcmctdem.exe garjmcmctdem.con
pause
