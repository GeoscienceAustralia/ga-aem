echo off

REM splitasciibycolumn.exe can be found in the open-source GitHub repository https://github.com/rcb547/utility-programs.
REM The Windows executables are available https://github.com/rcb547/utility-programs/releases/tag/v1.0.

CD output
splitasciibycolumn.exe inversion.output.dat 5
MKDIR lines
MOVE *.asc lines

pause


