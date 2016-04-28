GA-AEM Source Code Repository
=============================

Geoscience Australia Airborne Electromagnetics Programs

Author: Ross C Brodie, Geoscience Australia
	ross.c.brodie at ga.gov.au

Language:	C++

Third party software included in this repository
1.	TEMPLATE NUMERICAL TOOLKIT (TNT)
	http://math.nist.gov/tnt/index.html

2.	FFTW
	http://www.fftw.org/

Third party software dependencies
1.	MPI
2.	FFTW

Currently included programs

1.	GAFORWARDMODELTDEM - 1D forward modelling program for AEM data
	To build on linux
		cd makefiles
		make -f gaforwardmodeltdem.make allclean

2.	GALEISBSTDEM - deterministic 1D sample by sample inversion of AEM data
	To build on linux
		cd makefiles
		make -f galeisbstdem.make allclean

3.	GARJMCMCTDEM - stochastic 1D sample by sample inversion of AEM data
	To build on linux
		cd makefiles
		make -f garjmcmctdem.make allclean

