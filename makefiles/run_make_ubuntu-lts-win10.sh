#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export csv_include='../submodules/csv-parser/single_include'
export ticpp_dir='../submodules/ticpp'
export FFTW_DIR='/usr/lib/x86_64-linux-gnu'
export eigen_include='../submodules/eigen'
export petsc_include='/usr/include/petsc'

if [ "$compiler" = 'intel' ] ; then
	echo 'Building with Intel compiler'
	#module load intel-compiler
	export cxx=icpc
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas -Wno-unused-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security'
	export exedir='../bin/intel'
elif [ "$compiler" = 'gnu' ] ; then
	echo 'Building with GCC compiler'
	#module load gcc/system
	export cxx=g++
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas -Wno-unused-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security'
	export exedir='../bin/gnu'
else 
	echo 'Unknown compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_GDAL=1

if [ $HAVE_NETCDF = 1 ] ; then
	echo 'Building with NETCDF'
	export geophysics_netcdf_include='../submodules/geophysics-netcdf/src'
	export marray_include='../submodules/geophysics-netcdf/submodules/marray/include/andres'
fi

if [ $HAVE_GDAL = 1 ] ; then
	echo 'Building with GDAL'
	export gdal_include='/usr/include/gdal'
fi

#module list
echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme
echo HAVE_NETCDF = $HAVE_NETCDF
echo HAVE_GDAL = $HAVE_GDAL

echo ---------------------------------------

#Compiled as shared libs
#make -f gatdaem1d_python.make $makemode
#make -f gatdaem1d_matlab.make $makemode

#Compile without MPI
#make -f ctlinedata2sgrid.make $makemode
#make -f ctlinedata2slicegrids.make $makemode
#make -f example_forward_model.make $makemode
make -f gaforwardmodeltdem.make $makemode

#Compile with MPI
make -f galeisbstdem.make $makemode
#make -f garjmcmctdem.make $makemode
#make -f galeiallatonce.make $makemode
#make -f galeisbsfdem.make $makemode

