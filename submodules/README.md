Submodules
===========

Git submodules are used in this project. They contain utility and third party source code and binaries. Depending on what you are doing, you may not need all the submodules. When the repository is initially cloned the submodule directories will not be populated with code/files. To populate them you need to issue the "git submodule init" and "git submodule update" commands to initialise and populate the submodules with the (correct commit/version) of each repository. For example,

```bash script
>> cd myrepos/ga-aem 
>> git submodule init 
>> git submodule update 
>> git submodule status 
```

# cpp-utils
	- C++ utility code, classes and functions used in a variety of separate repositories.
	- Author: Ross C Brodie
	- Language: C++
	- Repository: git@github.com:rcb547/cpp-utils.git
	- Only required if you want to compile the source code.
	- Not required if you just want to run the precompiled Windows executables.

# tnt 
	- Template Numerical Toolkit (TNT): Linear Algebra Module
	- Authors: National Institute of Standards and Technology (NIST) by employees of the Federal Government in the course of their official duties.
	- Language: C++
	- Website: http://math.nist.gov/tnt/
	- Repository: git@github.com:rcb547/tnt.git
	- Not required if you just want to run the precompiled Windows executables.
	- Only required if you want to compile the source code.

# fftw3.2.2.dlls
	- Precompiled FFTW 3.2.2 Windows DLLs  
	- Downloaded from: http://www.fftw.org/install/windows.html  
	- Website: http://www.fftw.org   
	- Repository: git@github.com:rcb547/fftw3.2.2.dlls.git
	- Applies to Microsoft Windows only
	- Required if you want to execute the forward modelling or inversion programs on Windows.
	- Required if you want to compile the forward modelling or inversion programs on Windows.
