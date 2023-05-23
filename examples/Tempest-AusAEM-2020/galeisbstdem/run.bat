@echo off

REM Add executable and FFTW directories to your search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

REM Run standalone on a single CPU
galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con			&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets
REM galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con	&:: Invert X and Z-components total-field data and solve Dx and Dz Tx-Rx offsets and the Rx pitch
REM galeisbstdem.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con	&:: Invert only Z-component secondary field and do not solve for ant geometry

REM Run using 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con 4		&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets
REM galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con 4	&:: Invert X and Z-components total-field data and solve Dx and Dz Tx-Rx offsets and the Rx pitch
REM galeisbstdem.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con 4	&:: Invert only Z-component secondary field and do not solve for ant geometry

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-xzamplitude-solve-offsets.con			&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con	&:: Invert X and Z-components total-field data and solve Dx and Dz Tx-Rx offsets and the Rx pitch
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con	&:: Invert only Z-component secondary field and do not solve for ant geometry

REM If you do not have MPI installed substitute galeisbstdem-nompi.exe in place of galeisbstdem.exe
REM galeisbstdem-nompi.exe galeisbstdem-xzamplitude-solve-offsets.con			&:: Invert amplitude of total-field data in the X&Z-component plane and do not solve Dx and Dz Tx-Rx offsets
REM galeisbstdem-nompi.exe galeisbstdem-xzcomponents-solve-offsets-and-rxpitch.con	&:: Invert X and Z-components total-field data and solve Dx and Dz Tx-Rx offsets and the Rx pitch
REM galeisbstdem-nompi.exe galeisbstdem-zcomponentonly-do-not-solve-geometry.con	&:: Invert only Z-component secondary field and do not solve for ant geometry

pause


