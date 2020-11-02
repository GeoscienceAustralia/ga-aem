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
	module load intel-cc/12.1.9.293
	export cxx=icpc
	export cxxflags='-std=c++11 -O3 -Wall'
	export exedir='../bin/vdi/intel'
elif [ $compiler == 'gnu' ] ; then
	echo 'Building with GCC compiler'
	#module load gcc/9.1.0
	export cxx=g++
	export cxxflags='-std=c++11 -O3 -Wall -Wno-unknown-pragmas'
	export exedir='../bin/vdi/gnu'
else 
	echo 'Unknow compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_NETCDF=1
export HAVE_CGAL=0
module load openmpi/4.0.2
module load fftw3/3.3.8
#module load eigen/3.3.4
#module load petsc/3.9.4

if [ $HAVE_NETCDF == 1 ] ; then
	echo 'Building with NETCDF'
	module load netcdf/4.7.4
	#replace with your netCDF installation root.
	#not necessary if netCDF is already in your paths
	#export CPATH=/g/data/r78/rlt118/netCDF/include/:$CPATH
	#export LIBRARY_PATH=/g/data/r78/rlt118/netCDF/lib/:$LIBRARY_PATH
	#export LD_LIBRARY_PATH=/g/data/r78/rlt118/netCDF/lib/:$LD_LIBRARY_PATH
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
#make -f gatdaem1d_python.make $makemode
#make -f gatdaem1d_matlab.make $makemode

#Compile without MPI
#make -f ctlinedata2sgrid.make $makemode
#make -f ctlinedata2slicegrids.make $makemode
#make -f example_forward_model.make $makemode
#make -f gaforwardmodeltdem.make $makemode

#Compile with MPI
#make -f galeisbstdem.make $makemode
make -f garjmcmctdem.make $makemode
#make -f galeiallatonce.make $makemode
#make -f galeisbsfdem.make $makemode

