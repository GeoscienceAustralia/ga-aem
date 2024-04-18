# GA-AEM Python forward modelling interface

## Description
The Python (>=v3.5) interface consists of a C/C++ shared libarary (.so on Linux or .dll on Windows) called gatdaem1d which contains time-domain forward modelling and derivative functions which are called by the Python interpreter.

## Compiling and installing the C/C++ shared libraries
First the shared library needs to be built with CMake.  See one of the CMake build scripts in the root directory of the ga-aem source code repository. If you are only interested in the Python interface, you need only build the ***`python_bindings`*** target.

## Install directory contents
After being built successfully the install directory should contain,
- [ga-aem-install-dir]/python contains the package set up or installation function `setup.py`.
- [ga-aem-install-dir]/python/gatdaem1d contains the file `__init__.py` which is the package's Python classes and function code. It is also where the compiled shared library will reside after compilation.
	- On Linux the shared library is [ga-aem-install-dir]/python/gatdaem1d/gatdaem1d.so.
	- On Windows the shared library is [ga-aem-install-dir]/python/gatdaem1d/gatdaem1d.dll.
- [ga-aem-install-dir]/python/examples contains example Python usage code.

## PIP install of the Python package
- To install as a python package you can then,
```bash
	cd [ga-aem-install-dir]/python
	python -m pip install .
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

