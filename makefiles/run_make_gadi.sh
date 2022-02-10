#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export geophysics_netcdf_include='../submodules/geophysics-netcdf/src'
export marray_include='../submodules/geophysics-netcdf/submodules/marray/include/andres'
export csv_include='../submodules/csv-parser/single_include'

if [ $compiler == 'intel' ] ; then
	echo 'Building with Intel compiler'
	module load intel-compiler
	export cxx=icpc
	export cc=icc
	export ccflags='-O3 -Wall'
	export cxxflags='-std=c++17 -O3 -Wall'
	export exedir='../bin/gadi/intel'
elif [ $compiler == 'gnu' ] ; then
	echo 'Building with GCC compiler'
	module load gcc/11.1.0
	export cxx=g++
	export cc=gcc
	export ccflags='-O3 -Wall -Wno-unknown-pragmas'
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas'
	export exedir='../bin/gadi/gnu'
else 
	echo 'Unknow compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_GDAL=1
export HAVE_CGAL=0
module load openmpi/4.0.1
module load fftw3/3.3.8
module load eigen/3.3.7
module load petsc/3.12.2

if [ $HAVE_NETCDF == 1 ] ; then
	echo 'Building with NETCDF'
	export geophysics_netcdf_root='/home/547/rcb547/code/repos/geophysics-netcdf'
	module load netcdf/4.7.1
fi

if [ $HAVE_GDAL == 1 ] ; then
	echo 'Building with GDAL'
	module load gdal/3.0.2
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

#Compiled as shared libraries
#make -f gatdaem1d_python.make $makemode
#make -f gatdaem1d_matlab.make $makemode
#make -f gatdaem1d_julia.make  $makemode

#Compiled as static C-callable library
make -f gatdaem1d_c_library.make $makemode
#Compiled with C to use the C-callable library
make -f example_forward_model_c.make $makemode

#Compiled without MPI
#make -f ctlinedata2sgrid.make $makemode
#make -f ctlinedata2slicegrids.make $makemode
#make -f removelog10conductivityfromsgrid.make $makemode
#make -f example_forward_model.make $makemode
#make -f gaforwardmodeltdem.make $makemode

#Compiled with MPI
make -f galeisbstdem.make $makemode
make -f garjmcmctdem.make $makemode
#make -f galeiallatonce.make $makemode
#make -f galeisbsfdem.make $makemode

