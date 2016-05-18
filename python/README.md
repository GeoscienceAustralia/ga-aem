#ga-aem/python
======

  * Directory for the PYTHON(3) interface module gatdaem1d which contains time domain forward modelling and derivative functions.

##Module
  * ga-aem/python/gatdaem1d/
  * Contains the C/C++ code compiled into a shared library module  
  * The shared library (.so on linux or .dll on Windows) is called by PYTHON  
  * The ,pyc file(s) contain the Python classes and functions  

##Examples
  * ga-aem/python/examples
  * Contains examples of how to use the gatdaem1d module  
  * You will probably need to set you PYTHONPATH variable  
  * For example on Windows "set PYTHONPATH=Z:\code\repos\ga-aem\python;%PYTHONPATH%"

##Python dependencies
  * The gatdaem1d examples require the numpy and matplotlib packages.
  * Numpy and matplotlib can usually be installed as follows:
  * >> python -m pip install numpy
  * >> python -m pip install matplotlib

##Compiling the C/C++ shared libraries
  * Linux: On linux the shared library (.so) will need to be compiled 
  * see ga-aem/makefiles/gatdaem1d_python.make
  * Windows: On Windows a 64 bit Winodws dll is included and may work if you are using a 64 bit python interpreter
  * otherwise a Visual Studio 2013 project is located in ga-aem\vs2013\gatdaem1d_python

##FFTW
  * Before running the examples make sure the FFTW libraries are in your search path

