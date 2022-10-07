ga-aem/makefiles
================

##Directory for makefiles
	*To compile on linux: Read, modify and run version of "my_run_make.sh"
	*See for example raijin_run_make.sh

##You will need to set the
	* cxx (C++ compiler)
	* mpicxx (MPI compiler script)
	* cxxflags (C++ compiler flags),
	* and exedir (executable directory)

##For example:
	* #GNU GCC compiler on raijin.nci.org.au
	* export cxx=g++
	* export mpicxx=mpiCC
	* export cxxflags='-std=c++11 -O3 -Wall'
	* export exedir='../bin/raijin/gnu'


