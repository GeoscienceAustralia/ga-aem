# Submodules for ga-aem

# Description
Git submodules are used in the ga-aem git project. They contain utility and third party source code. Depending on what you are doing, you may not need all the submodules. They are not required if you just want to run the precompiled Windows executables.

# Initialisation and updating
By default, when the ga-aem repository is initially cloned the submodule directories will not be populated with code/files. To populate them, you need to issue the `git submodule init` and `git submodule update --recursive` commands to initialise and populate each submodule with the (correct commit/version) of each repository.  The `--recursive` switch is required because some submodules have their own submodules.  For example,
```bash
> cd <ga-aem-repository-directory>
> git submodule init
> git submodule update --recursive
> git submodule status 
```
Better still, when initially cloning the ga-aem repo you should use `--recursive` as follows to initialise and populate all submodules and their respective submodules.
```bash
> git clone --recursive https://github.com/GeoscienceAustralia/ga-aem.git
```

# The submodules

## cpp-utils
	- cpp-utils is a header-only source code library of classes and functions used in a variety of separate repositories.
	- Author: Ross C Brodie
	- Language: C++
	- Repository: git@github.com:rcb547/cpp-utils.git
	- Not required if you just want to run the precompiled Windows executables.
## geophysics-netcdf
	- geophysics-netcdf is a header source code library that implements reading and writing Geoscience Australia's NETCDF geophysical line data format.
	- Author: Ross C Brodie
	- Language: C++
	- Repository: git@github.com:GeoscienceAustralia/geophysics-netcdf.git
	- Not required if you just want to run the precompiled Windows executables.
## eigen
	- Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
	- Language: C++
	- Website: http://eigen.tuxfamily.org/.	
	- Repository: https://github.com/eigenteam/eigen-git-mirror.git	
	- Not required if you just want to run the precompiled Windows executables.
	- Only required if you want to compile the source code.
## tinyxml
	- TiCPP is short for the official name TinyXML++ and is a simple XML reader/writer.
	- Authors: 	Ryan Pusztai, Ryan Mulder
	- Language: C++
	- Website: https://github.com/wxFormBuilder/ticpp
	- Repository: git@github.com:wxFormBuilder/ticpp.git
	- Not required if you just want to run the precompiled Windows executables.
	- Only required if you want to compile the source code for ctlinedatatocurtainimage.
## netcdf-cxx4
	- netcdf-cxx4 is a source code library that implements the C++ binding for the NETCDF C library.
	- Author: Lynton Appel
	- Language: C++
	- Website: https://github.com/Unidata/netCDF-cxx4
	- Repository: git@github.com:Unidata/netcdf-cxx4.git
## csv-parser
	- CSV-Parser is a CSV file parsing library.
	- Author: Vincent La
	- Language: C++
	- Website: https://github.com/vincentlaucsb/csv-parser
	- Repository: https://github.com/vincentlaucsb/csv-parser.git
	- included as a submodule of the cpp-utils submodule
## marray
	- Fast Runtime-Flexible Multi-dimensional Arrays and Views in C++
	- Author: Bjoern Andres
	- Language: C++
	- Website: http://www.andres.sc/marray.html
	- Repository: https://github.com/bjoern-andres/marray.git
	- included as a submodule of the geophysics-netcdf submodule
