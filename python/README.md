# GA-AEM Python forward modelling interface

## Description
The Python (>=v3.5) interface consists of a C/C++ shared library (.so on Linux or .dll on Windows) called gatdaem1d which contains time-domain forward modelling and derivative functions which are called by the Python interpreter.

## Compiling and installing the C/C++ shared libraries
Ensure that the GNU C++ compiler is available, and FFTW has been installed.

Mac OSX currently has issues using cmake and brew installed gcc compilers. So we are stuck using a Makefile. Simply type "make" in the python folder.  On Linux, the easiest option is to use "make" also.

Two environment variables need to be set before compilation using "export cxx=g++" and "export FFTW_ROOT=<path to FFTW>".

On Windows, follow the documentation in the root folder.  Once the library is compiled, go ahead and pip install.

## PIP install of the Python package
- To install as a python package you can then,
```bash
	pip install .
```
- Note that **`python`** may need to be **`python3`** on your system, and make sure the correct environment is activated/sourced.

## Examples
- The directory [ga-aem-install-dir]/python/examples contains an example of how to use the gatdaem1d package.
- To run the example
```bash
	cd [ga-aem-install-dir]/python/examples
	python skytem_example.py
```

## Dependencies

### Python dependencies
- The gatdaem1d package itself requires the ***numpy*** Python package.
- The example additionally requires the ***matplotlib*** Python package.
- Numpy and matplotlib can usually be installed as follows:
````bash
	python -m pip install numpy
	python -m pip install matplotlib
````

### FFTW dependency
- Before running the examples make sure the FFTW libraries are in your system's search path

