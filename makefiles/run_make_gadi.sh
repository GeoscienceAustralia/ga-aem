#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

module load openmpi/4.1.4
module load fftw3/3.3.8
module load eigen/3.3.7
module load petsc/3.17.4

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export csv_include='../submodules/csv-parser/single_include'
export ticpp_dir='../submodules/ticpp'
#export FFTW_DIR='/usr/lib/x86_64-linux-gnu'
export eigen_include='../submodules/eigen'
export petsc_include='/usr/include/petsc'


if [ $compiler == 'intel' ] ; then
	echo 'Building with Intel compiler'
	module load intel-compiler
	export cxx=icpc
	export cxxflags='-std=c++17 -O3  -Wall -Wno-unknown-pragmas -Wno-unused-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security'
	export exedir='../bin/gadi/intel'
elif [ $compiler == 'gnu' ] ; then
	echo 'Building with GCC compiler'	
	module load gcc/11.1.0
	export cxx=g++
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas -Wno-unused-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security'
	export exedir='../bin/gadi/gnu'
else 
	echo 'Unknow compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_GDAL=1

if [ $HAVE_NETCDF == 1 ] ; then
	echo 'Building with NETCDF'	
	module load netcdf/4.8.0
	export geophysics_netcdf_include='../submodules/geophysics-netcdf/src'
	export marray_include='../submodules/geophysics-netcdf/submodules/marray/include/andres'
fi

if [ $HAVE_GDAL == 1 ] ; then
	echo 'Building with GDAL'
	module load gdal/3.0.2
	export gdal_include='/usr/include/gdal'
fi

module list
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
#make -f gatdaem1d_julia.make  $makemode

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

