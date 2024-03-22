# Python
Directory containing the Python(3) interface files package.

## Description
The Python (>=v3.5) interface consists of a C/C++ shared libarary (.so on Linux or .dll on Windows) called gatdaem1d which contains time domain forward modelling and derivative functions which are called by the Python interpreter.
- [ga-aem-install-dir]/python contains the package set up or installation function setup.py.
- [ga-aem-install-dir]/python/gatdaem1d contains the file \_\_init\_\_.py which is the package's Python classes and function code.  It is also where the compiled shared library will reside after compilation.
- [ga-aem-install-dir]/python/examples/ contains example Python usage code.

## Compiling the C/C++ shared libraries
- First the shared library needs to be built with CMAKE.  See one of the installation scripts in the root directory of the repository.
- After being built successfully,
	- The shared library are written to [ga-aem-install-dir]/python/gatdaem1d.
	- On Linux the shared library is [ga-aem-install-dir]/python/gatdaem1d/gatdaem1d.so.
	- On Windows the shared library is [ga-aem-install-dir]/python/gatdaem1d/gatdaem1d.dll.
	- On Windows a 64 bit Winodws dll is included in the repository and may work out of the box if you are using a 64 bit Python interpreter.

## Installation of the Python package
- To install the package you can then,
```bash
	cd [ga-aem-install-dir]/python
	python -m pip install .
```

## Examples
- The directory [ga-aem-install-dir]/python/examples contains an example of how to use the gatdaem1d package.
- To run the example
```bash
	cd [ga-aem-install-dir]/examples
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

