Submodules
===========

Git submodules are used in this project. They contain utility and third party source code and binaries. Depending on what you are doing, you may not need all the submodules.

# Initialisation and updating

By default, when the repository is initially cloned the submodule directories will not be populated with code/files. To populate them you need to issue the "git submodule init" and "git submodule update --recursive" commands to initialise and populate the submodules with the (correct commit/version) of each repository.  The "--recursive" is required because some submodules have their own submodules.  For example,

```script
>> cd myrepos/ga-aem 
>> git submodule init
>> git submodule update --recursive
>> git submodule status 
```

Alternatively, when initially cloning the repo you can use "--recursive" as follows to initialise and populate all submodules and their respective submodules.

```script
>> git clone --recursive git@github.com:rcb547/ga-aem-csiro.git
```


# The submodules

## cpp-utils
	- C++ utility code, classes and functions used in a variety of separate repositories.
	- Author: Ross C Brodie
	- Language: C++
	- Repository: git@github.com:rcb547/cpp-utils.git
	- Only required if you want to compile the source code.
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
	- Language: C++
	- Website: https://github.com/wxFormBuilder/ticpp
	- Repository: git@github.com:wxFormBuilder/ticpp.git
	- Not required if you just want to run the precompiled Windows executables.
	- Only required if you want to compile the source code for ctlinedatatocurtainimage.

## csv-parser
	- CSV-Parser is a CSV file parsing library.
	- Language: C++
	- Website: https://github.com/vincentlaucsb/csv-parser
	- Repository: https://github.com/vincentlaucsb/csv-parser.git

## geophysics-netcdf
	- geophysics-netcdf is a header source code library that implements reading and writing Geoscience Australia's NETCDF geophysical line data format.
	- Author: Ross C Brodie
	- Language: C++
	- Repository: git@github.com:GeoscienceAustralia/geophysics-netcdf.git

## netcdf-cxx4
	- netcdf-cxx4 is a source code library that implements the C++ binding for the NETCDF C library.
	- Author: Lynton Appel
	- Language: C++
	- Website: https://github.com/Unidata/netCDF-cxx4
	- Repository: git@github.com:Unidata/netcdf-cxx4.git
