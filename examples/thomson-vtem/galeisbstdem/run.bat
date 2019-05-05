@echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)

REM set path=..\..\..\..\fftw3.2.2.dlls\64bit;%path%
set path=C:\Microsoft_HPC_Pack_2012\Bin;%path%
set path=..\..\..\bin\x64\Release\;%path%

REM Note: The control files galeisbstdem_ascii_dfn.con and galeisbstdem_ascii_column.con have the same inversion settings. However the _dfn.con version uses the ASEGGDF2 .dfn file variable names and the _column.con version uses column numbers to define fields.
REM e.g.  In galeisbstdem_ascii_column.con TX_Height = Column 30
REM       In galeisbstdem_ascii_dfn.con    TX_Height = emloop_height

REM Use a single standalone CPU
galeisbstdem.exe galeisbstdem_ascii_dfn.con
REM galeisbstdem.exe galeisbstdem_ascii_column.con

REM Use 4 MPI processes
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem_ascii_dfn.con
REM mpiexec -np 4 galeisbstdem.exe galeisbstdem_ascii_column.con

REM Use 4 OpenMP threads
REM galeisbstdem.exe galeisbstdem_ascii_dfn.con 4
REM galeisbstdem.exe galeisbstdem_ascii_column.con 4

REM If you don't have MPI installed use the non-MPI version
REM Use 4 OpenMP threads 
REM galeisbstdem-nompi.exe galeisbstdem_ascii_dfn.con 4
REM galeisbstdem-nompi.exe galeisbstdem_ascii_column.con 4


pause

