#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export ompmode=$3
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'

if [ $compiler = 'intel' ] ; then
	echo 'Building with Intel compiler'
	export cxx=icpc
	export cxxflags='-std=c++11 -O3 -Wall'
	export exedir='../bin/ubuntu/intel'
	if [ $ompmode = 'openmp' ] ; then
		echo 'Compiling with openMP multithreading'
		export cxxflags="${cxxflags} -openmp"
		export ldflags="-fopenmp"
	fi
elif [ $compiler = 'gnu' ] ; then
	echo 'Building with GCC compiler'
	export cxx=g++
	export cxxflags='-std=c++11 -O3 -Wall -Wno-unknown-pragmas'
	export exedir='../bin/ubuntu/gnu'
	if [ $ompmode = 'openmp' ] ; then
		echo 'Compiling with openMP multithreading'
		export cxxflags="${cxxflags} -fopenmp"
		export ldflags="-fopenmp"
	fi
elif [ $compiler = 'gnu-debug' ] ; then
	echo 'Building with GCC compiler with debug symbols'
	export cxx=g++
	export cxxflags='-std=c++11 -Og -Wall -Wno-unknown-pragmas -g'
	export exedir='../bin/ubuntu/gnu'
	if [ $ompmode = 'openmp' ] ; then
		echo 'Compiling with openMP multithreading'
		export cxxflags="${cxxflags} -fopenmp"
		export ldflags="-fopenmp"
	fi
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
