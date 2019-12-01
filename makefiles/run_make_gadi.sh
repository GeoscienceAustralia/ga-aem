#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'

if [ $compiler == 'intel' ] ; then
	echo 'Building with Intel compiler'
	module load intel-compiler
	export cxx=icpc
	export cxxflags='-std=c++11 -O3 -Wall'
	export exedir='../bin/gadi/intel'
elif [ $compiler == 'gnu' ] ; then
	echo 'Building with GCC compiler'
	module load gcc/system
	export cxx=g++
	export cxxflags='-std=c++11 -O3 -Wall -Wno-unknown-pragmas'
	export exedir='../bin/gadi/gnu'
else 
	echo 'Unknow compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_GDAL=0
export HAVE_CGAL=0
module load openmpi/4.0.1
module load fftw3/3.3.8
module load eigen/3.3.7

if [ $HAVE_NETCDF == 1 ] ; then
	echo 'Building with NETCDF'
	export geophysics_netcdf_root='/home/547/rcb547/code/repos/geophysics-netcdf'
	module load netcdf/4.7.1
fi

if [ $HAVE_GDAL == 1 ] ; then
	echo 'Building with GDAL'
	module load gdal
fi

if [ $HAVE_CGAL == 1 ] ; then
	echo 'Building with CGAL'
	module load cgal
fi


module list
echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme

echo ---------------------------------------

#Compiled as shared libs
make -f gatdaem1d_python.make $makemode
make -f gatdaem1d_matlab.make $makemode

#Compile without MPI
make -f ctlinedata2sgrid.make $makemode
# Need GDAL make -f ctlinedata2slicegrids.make $makemode
make -f example_forward_model.make $makemode
make -f gaforwardmodeltdem.make $makemode

#Compile with MPI
make -f galeisbstdem.make $makemode
make -f garjmcmctdem.make $makemode
# make -f galeisbsfdem.make $makemode
# NEED PETSC make -f galeiallatonce.make $makemode


