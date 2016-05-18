python
======

Directory for the PYTHON interface module gatdaem1d which contains time domain forward modelling and derivative functions.

/python/gatdaem1d/ contains the module for the time-domain shared library (it is a shared library .so or .dll that PYTHON can call)
/python/examples contains an example of how to use gatdaem1d module

Make sure the FFTW dlls (ga-aem\third_party\fftw3.2.2.dlls\64bit) are in your Windows path variable before running the examples.

In case the shared libraries need to be regenerated see:
    Windows:	A Visual Studio 2013 project is located in ga-aem\vs2013\gatdaem1d_python
    Linux:	ga-aem/makefiles/gatdaem1d_python.make

