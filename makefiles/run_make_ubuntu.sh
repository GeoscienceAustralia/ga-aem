#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export geophysics_netcdf_include='../submodules/geophysics-netcdf/src'
export eigen_include='../submodules/eigen'
export marray_include='../submodules/geophysics-netcdf/submodules/marray/include/andres'
export csv_include='../submodules/csv-parser/single_include'

export FFTW_DIR='/usr/lib/x86_64-linux-gnu'

if [ $compiler = 'intel' ] ; then
	echo 'Building with Intel compiler'
	export cxx=icpc
	export cc=icc
	export cxxflags='-std=c++17 -O3 -Wall'
	export exedir='../bin/ubuntu/intel'
elif [ $compiler = 'gnu' ] ; then
	echo 'Building with GCC compiler'
	export cxx=g++
	export cc=gcc
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas -Wno-unused-variable -Wno-sign-compare'
	export exedir='../bin/ubuntu/gnu'
else 
	echo 'Unknown compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_CGAL=0

if [ $HAVE_NETCDF = 1 ] ; then
	echo 'Building with NETCDF'
fi

if [ $HAVE_CGAL = 1 ] ; then
	echo 'Building with CGAL'
fi


echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme

echo ---------------------------------------

#Compiled as shared libs
#make -f gatdaem1d_python.make $makemode
#make -f gatdaem1d_matlab.make $makemode

#Compiled as static C-callable library
#make -f gatdaem1d_c_library.make $makemode
#Compiled with C to use the C-callable library
#make -f example_forward_model_c.make $makemode


#Compile without MPI
#make -f ctlinedata2sgrid.make $makemode
#make -f ctlinedata2slicegrids.make $makemode
#make -f example_forward_model.make $makemode
#make -f gaforwardmodeltdem.make $makemode

#Compile with MPI
make -f galeisbstdem.make $makemode
#make -f garjmcmctdem.make $makemode
#make -f galeiallatonce.make $makemode
#make -f galeisbsfdem.make $makemode

