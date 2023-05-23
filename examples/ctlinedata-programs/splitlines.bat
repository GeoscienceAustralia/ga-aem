echo off

CD output
splitasciibycolumn.exe inversion.output.dat 5
MKDIR lines
MOVE *.asc lines

pause


