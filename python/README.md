# GA-AEM Python forward modelling interface

## Description
The Python (>=v3.5) interface consists of a C/C++ shared library (.so on Linux or .dll on Windows) called gatdaem1d which contains time-domain forward modelling and derivative functions which are called by the Python interpreter.

## Compiling and installing the C/C++ shared libraries
We are using Meson to automatically build the required C++ library.
Ensure that the GNU C++ compiler is available, and FFTW has been installed.

## PIP install of the Python package
- To install as a python package you can then,
```bash
	cd [ga-aem-install-dir]/python
	pip install .
```
- Note that **`python`** may need to be **`python3`** on your system.

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

