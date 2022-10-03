Matlab
======

Directory for the MATLAB interface module (C callable shared libarary) gatdaem1d which contains time domain forward modelling and derivative functions.
	- matlab/bin/libgatdaem1d.so is the time-domain Linux shared library
	- matlab/bin/gatdaem1d.mexw64 is the time-domain 64-Bit Windows shared library (it is a dll that MATLAB can call)
	- matlab/gatdaem1d_functions contains the wrapper functions MATLAB scripts .m
	- matlab/rjmcmctdem_functions contains the wrapper functions MATLAB scripts .m for handling Monte-Carlo invesrion outputs
	- matlab/examples contains examples of how to use gatdaem1d module

Make sure the FFTW dlls (<ga-aem install dir>\third_party\fftw3.2.2.dlls\64bit) are in your Windows path variable before running the examples.

You should run 'mex -setup' from the Matlab command prompt first to allow Matlab to find your installed C Compiler.

In your Matlab environment add the following directories to your Matlab search path
	- ga-aem-install-dir/matlab/bin
	- ga-aem install dir/gatdaem1d_functions
	- ga-aem install dir/rjmcmctdem_functions

On Windows you may need to start Matlab from the "x64 Native Tools Command Prompt for VS 2019" (or similar Start Menu Item), or at least add the directory that contains the compiler executable "cl.exe" in your windows search path. For example add something like one of the following to your windows search path
	- C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64\cl.exe
	- C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\bin\Hostx64\x64\cl.exe
