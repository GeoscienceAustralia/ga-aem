echo off

REM Add executable and FFTW directories to yo search path
REM (Ideally you would add these to your PATH environment variable)
set path=..\..\..\bin\x64\Release\;%path%
set path=..\..\..\third_party\fftw3.2.2.dlls\64bit;%path%

gaforwardmodeltdem.exe skytem_lm.con
gaforwardmodeltdem.exe skytem_hm.con

pause
