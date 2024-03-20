# Matlab

Directory for the MATLAB interface module (C callable shared libarary) gatdaem1d which contains time domain forward modelling and derivative functions.

# Contents
- matlab/bin/libgatdaem1d.so is the time-domain Linux shared library
- matlab/bin/gatdaem1d.mexw64 is the time-domain 64-Bit Windows shared library (it is a dll that MATLAB can call)
- matlab/gatdaem1d_functions contains the wrapper functions MATLAB scripts .m
- matlab/rjmcmctdem_functions contains the wrapper functions MATLAB scripts .m for handling Monte-Carlo invesrion outputs
- matlab/examples contains examples of how to use the gatdaem1d module

# Notes
1. You should run 'mex -setup' from the Matlab command prompt first to allow Matlab to find your installed C Compiler.
2. On Windows,
	- you will need to have Microsoft Visual Studio 2019 or 2022 installed.
	- you will probably need to start\launch Matlab from the "x64 Native Tools Command Prompt for VS 2022" (or similar Start Menu Item).
	- make sure that the directory containing your installed FFTW libraries (e.g. %LocalAppData%\fftw-3.3.5-dll64) is in your PATH environment variable before running the examples.  That folder will contain libfftw3-3.dll.
3. In your Matlab environment add the following directories to your Matlab search path:
	- [ga-aem-install-dir]/matlab/bin
	- [ga-aem install dir]/gatdaem1d_functions
	- [ga-aem install dir]/rjmcmctdem_functions

# Examples
There are several example of how to use the Matlab interface in the directory [ga-aem-install-dir]\matlab\examples
- The examples have paths setup to be run as if your Matlab current working directory is [ga-aem-install-dir]\matlab\examples
- If they are not already added in your Matlab startup.m script, you may need to add the gatdaem1d wrapper function (.m) files and the shared library to you Matlab path at the beginning of each script, for example you may need to uncomment and alter these lines,
```matlab
	addpath('..\bin');
	addpath('..\gatdaem1d_functions');
	addpath('C:\fftw-3.3.5-dll64');
```