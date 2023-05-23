@echo off

REM Add executable and FFTW directories to your search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM Run standalone on a single CPU
galeisbstdem.exe galeisbstdem.con	&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets

REM Run using 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem.con 4	&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem.con	&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets

REM If you do not have MPI installed substitute galeisbstdem-nompi.exe in place of galeisbstdem.exe
REM galeisbstdem-nompi.exe galeisbstdem.con		&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets

pause


