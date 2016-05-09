matlab
======

Directory for the MATLAB interface module gatdaem1d which contains time domain forward modelling and derivative functions.

matlab/bin/x64/gatdaem1d.mexw64 is the time-domain 64 Bit Windows shared library (it is a dll that MATLAB can call)
matlab/gatdaem1d_functions contains the wrapper functions MATLAB scripts .m
matlab/examples contains examples of how to use gatdaem1d module

Make sure the FFTW dlls (ga-aem\third_party\fftw3.2.2.dlls\64bit) are in your Windows path variable before running the examples.
A Visual Studio project is located in ga-aem\vs2013\gatdaem1d_matlab in case the mexw64 needs to be recreated.

