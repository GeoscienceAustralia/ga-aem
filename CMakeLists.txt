cmake_minimum_required(VERSION 3.12)

# Set project name
project(ga-aem)

message(STATUS Build Type = ${CMAKE_BUILD_TYPE})

# Set options
option(WITH_MPI "WITH_MPI" ON)
option(WITH_NETCDF "WITH_NETCDF" ON)
option(WITH_GDAL "WITH_GDAL" ON)
option(WITH_PETSC "WITH_PETSC" ON)

message(STATUS WITH_MPI = ${WITH_MPI})
message(STATUS WITH_NETCDF = ${WITH_NETCDF})
message(STATUS WITH_GDAL = ${WITH_GDAL})
message(STATUS WITH_PETSC = ${WITH_PETSC})

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Set the compile flags
if("${CMAKE_C_COMPILER_ID}" MATCHES "MSVC" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   /O2")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")

	# Adding filesystem library
	link_libraries(-lstdc++fs)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -Wno-date-time")
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-unused-result -Wno-date-time")
	#On GCC, even with -Wno-date-time, still get warings of the form: warning: macro "__DATE__" might prevent reproducible builds [-Wdate-time]
endif()



# Configure FFTW
#find_package(FFTW REQUIRED fftw3 IMPORTED_TARGET)
#include_directories(PkgConfig::FFTW)
#link_libraries     (PkgConfig::FFTW)
#if(FFTW_FOUND)
#	message(STATUS "FFTW ${FFTW_VERSION} was found")
#endif()

#$find_package(FFTW IMPORTED_TARGET)
#if(FFTW_FOUND)
#	message(STATUS "FFTW ${FFTW_VERSION} was found")
#else()
#	set (FFTW_INCLUDE_DIRS ${FFTW_ROOT})
#	set (FFTW_LIBRARIES ${FFTW_ROOT})
#endif()

set(FFTW_ROOT $ENV{FFTW_ROOT})
#set(FFTW_INCLUDE_DIR $ENV{FFTW_INCLUDE_DIR})
#set(FFTW_LIBRARY_DIR $ENV{FFTW_LIBRARY_DIR})
#message(STATUS "FFTW include dir: ${FFTW_INCLUDE_DIR}")
#message(STATUS "FFTW library dir: ${FFTW_LIBRARY_DIR}")

find_path(FFTW_INCLUDE_DIRS
        NAMES "fftw3.h"
        PATHS ${FFTW_ROOT}
        PATH_SUFFIXES "include"
        NO_DEFAULT_PATH)

find_library(
        FFTW_LIB
        NAMES "fftw3-3" libfftw3-3
        PATHS ${FFTW_ROOT}
        PATH_SUFFIXES "lib" "lib64"
        NO_DEFAULT_PATH)


message(STATUS "FFTW include dirs: ${FFTW_INCLUDE_DIRS}")
message(STATUS "FFTW library: ${FFTW_LIB}")

# Configure MPI if opted for
if(${WITH_MPI})
	find_package(MPI)
	if(MPI_FOUND)
		message(STATUS "MPI was found")
		#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
		#add_definitions(-D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
		#link_libraries(${MPI_C_LIBRARIES})
	endif()
endif()

# Add OpenMP support (This is NOT the same as MPI)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND) 
	message(STATUS "OpenMP was found")
endif()

# Configure GDAL
if(${WITH_GDAL})
	find_package(GDAL REQUIRED)
	if(GDAL_FOUND)
		message(STATUS "GDAL ${GDAL_VERSION} was found")
		message(STATUS "GDAL_INCLUDE_DIRS ${GDAL_INCLUDE_DIRS}")
		message(STATUS "GDAL_LIBRARIES ${GDAL_LIBRARIES}")
	endif()
endif()

# Configure NETCDFCXX if opted for
if(${WITH_NETCDF_OLD})
	find_package(NETCDF REQUIRED) 
	if(NETCDF_FOUND)
		message(STATUS "NETCDF was found")
		message(STATUS "${netCDF_INSTALL_PREFIX}")
		message(STATUS "NETCDF install dir: ${netCDF_INSTALL_PREFIX}")
		message(STATUS "NETCDF include dir: ${netCDF_INCLUDE_DIR}")
		message(STATUS "NETCDF lib dir: ${netCDF_LIB_DIR}")
		set(netCDF_LIBRARIES netcdf hdf5 hdf5_hl curl_imp)
		#set(netCDF_LIBRARIES netcdf)
		message(STATUS "NETCDF libraries: ${netCDF_LIBRARIES}")
		#add_definitions(-DHAVE_NETCDF)
	endif()
endif()

if(${WITH_NETCDF})
	set(netCDF_INSTALL_PREFIX "C:/Program Files/netCDF 4.6.1")
	set(netCDF_INCLUDE_DIR ${netCDF_INSTALL_PREFIX}/include)
	set(netCDF_LIB_DIR ${netCDF_INSTALL_PREFIX}/lib)
	set(netCDF_LIBRARIES 
		${netCDF_LIB_DIR}/netcdf.lib 
		${netCDF_LIB_DIR}/hdf5.lib 
		${netCDF_LIB_DIR}/hdf5_hl.lib 
		${netCDF_LIB_DIR}/libcurl_imp.lib)

	#add_definitions(-DHAVE_NETCDF)

	message(STATUS "NETCDF install dir: ${netCDF_INSTALL_PREFIX}")
	message(STATUS "NETCDF include dir: ${netCDF_INCLUDE_DIR}")
	message(STATUS "NETCDF lib dir: ${netCDF_LIB_DIR}")
	message(STATUS "NETCDF libraries: ${netCDF_LIBRARIES}")
endif()

# Configure PETSc (only used for galeiallatonce)
if(${WITH_PETSC})	
	find_package(PETSC)
	if(PETSC_FOUND)		
		message(STATUS "PETSc was found")
		if (PETSC_INCLUDE_DIRS STREQUAL "")
			# If PETSc is not set up properly try a manual approach
			message(STATUS "PETSC pkg-config not set up properly -- resorting to manual setup of PETSC_INCLUDE_DIRS and PETSC_LDFLAGS")
			set(PETSC_INCLUDE_DIRS /apps/petsc/3.12.2/include)		
			set(PETSC_LDFLAGS "-L/apps/petsc/3.12.2/lib/ompi3/Intel -lpetsc -L/apps/intel-ct/2019.5.281/mkl/lib/intel64 -L/apps/intel-ct/2019.5.281/mkl/../compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/8 -L/apps/netcdf/4.7.1p/lib -L/apps/hdf5/1.10.5p/lib -L/apps/fftw3/3.3.8/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lsuperlu -lsuperlu_dist -lml -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lirc -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lnetcdf -lparmetis -lmetis -ltriangle -lm -lz -lX11 -lctetgen -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lrt -lquadmath -lstdc++ -ldl")
			message(STATUS "    PETSC_INCLUDE_DIRS=${PETSC_INCLUDE_DIRS}")	
			message(STATUS "    PETSC_LDFLAGS=${PETSC_LDFLAGS}")	
		endif()		
	else()
		message(STATUS "PETSc was NOT found -- will not build ${target}")
	endif()
endif()

if(${WITH_PETSC_OLD})	
	pkg_search_module(PETSC PETSc IMPORTED_TARGET) 
	if(PETSC_FOUND)		
		message(STATUS "PETSc was found")
		if (PETSC_INCLUDE_DIRS STREQUAL "")
			# If PETSc is not set up properly try a manual approach
			message(STATUS "PETSC pkg-config not set up properly -- resorting to manual setup of PETSC_INCLUDE_DIRS and PETSC_LDFLAGS")
			set(PETSC_INCLUDE_DIRS /apps/petsc/3.12.2/include)		
			set(PETSC_LDFLAGS "-L/apps/petsc/3.12.2/lib/ompi3/Intel -lpetsc -L/apps/intel-ct/2019.5.281/mkl/lib/intel64 -L/apps/intel-ct/2019.5.281/mkl/../compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/8 -L/apps/netcdf/4.7.1p/lib -L/apps/hdf5/1.10.5p/lib -L/apps/fftw3/3.3.8/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lsuperlu -lsuperlu_dist -lml -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lirc -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lnetcdf -lparmetis -lmetis -ltriangle -lm -lz -lX11 -lctetgen -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lrt -lquadmath -lstdc++ -ldl")
			message(STATUS "    PETSC_INCLUDE_DIRS=${PETSC_INCLUDE_DIRS}")	
			message(STATUS "    PETSC_LDFLAGS=${PETSC_LDFLAGS}")	
		endif()		
	else()
		message(STATUS "PETSc was NOT found -- will not build ${target}")
	endif()
endif()




# Add the include directories
include_directories(src/)
include_directories(submodules/cpp-utils/src/)
include_directories(submodules/csv-parser/single_include/)
include_directories(submodules/eigen/)
if(${WITH_NETCDF})
	include_directories(libnetcdfcxx INTERFACE ${netCDF_INCLUDE_DIR})
	include_directories(submodules/geophysics-netcdf/src/)
	include_directories(submodules/geophysics-netcdf/submodules/marray/include/andres/)
endif()


### Add the targets

# Add the NetCDF CXX library
file(GLOB NC_CXX_HEADERS submodules/netcdf-cxx4/cxx4/*.h)
file(GLOB NC_CXX_SOURCES submodules/netcdf-cxx4/cxx4/nc*.cpp)
add_library(libnetcdfcxx STATIC ${NC_CXX_SOURCES})
target_include_directories(libnetcdfcxx PUBLIC submodules/netcdf-cxx4/cxx4)
target_include_directories(libnetcdfcxx PUBLIC ${netCDF_INCLUDE_DIR})
target_link_directories(libnetcdfcxx PUBLIC ${netCDF_LIB_DIR})
target_link_libraries(libnetcdfcxx PUBLIC ${netCDF_LIBRARIES}) 
set(NETCDFCXX_INCLUDE_DIRS submodules/netcdf-cxx4/cxx4 ${netCDF_INCLUDE_DIR})
message(STATUS "NETCDFCXX_INCLUDE_DIRS ${NETCDFCXX_INCLUDE_DIRS}")

# Add the TICPP library submodule
add_subdirectory(submodules/ticpp)




#################
#################
#################

# Add the gatdaem1d static library
set(target gatdaem1d-static)
add_library(${target} STATIC src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
target_sources(${target} PRIVATE submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LIB})
install(TARGETS ${target} DESTINATION lib)

# Add the gatdaem1d shared library
set(target gatdaem1d-shared)
add_library(${target} SHARED src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
target_sources(${target} PRIVATE submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LIB})
#Install the libs
install(TARGETS ${target} DESTINATION bin)
install(TARGETS ${target} DESTINATION python)
install(TARGETS ${target} DESTINATION matlab)
#Install the include files
install(FILES   src/gatdaem1d.h TYPE INCLUDE)

# Add gaforwardmodeltdem executable
set(target gaforwardmodeltdem)
add_executable(${target} src/gaforwardmodeltdem.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LIB})
install(TARGETS ${target} DESTINATION bin)

# Add example_forward_model executable
set(target example_forward_model)
add_executable(${target} src/example_forward_model.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LIB})
#target_link_libraries(${target} PRIVATE cpp-utils-static)
install(TARGETS ${target} DESTINATION bin)


# Add the pure C example executable
set(target example_forward_model_c)
add_executable(${target} src/example_forward_model_c.c)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE gatdaem1d-static)
install(TARGETS ${target} DESTINATION bin)

# Add galeisbstdem executable
set(target galeisbstdem)
add_executable(${target} src/galeisbstdem.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LIB})
if(${WITH_MPI})
	target_compile_definitions(${target} PRIVATE -D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
	target_include_directories(${target} PRIVATE SYSTEM ${MPI_INCLUDE_PATH})
	target_link_libraries(${target} PRIVATE ${MPI_C_LIBRARIES})
endif()
if(${WITH_NETCDF})
	target_include_directories(${target} PRIVATE ${NETCDFCXX_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE libnetcdfcxx)
endif()
if(OpenMP_CXX_FOUND) 
	target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)	
endif()
install(TARGETS ${target} DESTINATION bin)

# Add garjmcmctdem executable
set(target garjmcmctdem)
if(${WITH_NETCDF})
	add_executable(${target} src/garjmcmctdem.cpp)
	target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
	target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${FFTW_LIB})
	target_include_directories(${target} PRIVATE ${NETCDFCXX_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE libnetcdfcxx)
	if(${WITH_MPI})
		target_compile_definitions(${target} PRIVATE -D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
		target_include_directories(${target} PRIVATE SYSTEM ${MPI_INCLUDE_PATH})
		target_link_libraries(${target} PRIVATE ${MPI_C_LIBRARIES})
	endif()
	install(TARGETS ${target} DESTINATION bin)
else()
	message(WARNING "${target} requires NETCDF - will not be built")
endif()

# Add galeiallatonce executable
set(target galeiallatonce)
if(MPI_FOUND AND PETSC_FOUND) 
	add_executable(${target} src/galeiallatonce.cpp)
	target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
	target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${FFTW_LIB})	
	target_include_directories(${target} PRIVATE ${PETSC_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${PETSC_LDFLAGS})
	install(TARGETS ${target} DESTINATION bin)
else()
	message(WARNING "${target} requires PETSc and MPI -- will not be built")
endif()

# Add galeisbsfdem executable
set(target galeisbsfdem)
add_executable(${target} src/galeisbsfdem.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
if(${WITH_MPI})
		target_compile_definitions(${target} PRIVATE -D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
		target_include_directories(${target} PRIVATE SYSTEM ${MPI_INCLUDE_PATH})
		target_link_libraries(${target} PRIVATE ${MPI_C_LIBRARIES})
endif()
if(OpenMP_CXX_FOUND) 
	target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)	
endif()
install(TARGETS ${target} DESTINATION bin)

# Add ctlinedata2georefimage executable
set(target ctlinedata2georefimage)
add_executable(${target} src/ctlinedata2georefimage.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
install(TARGETS ${target} DESTINATION bin)

# Add ctlinedata2curtainimage executable
if(${WITH_GDAL} AND GDAL_FOUND) 
	set(target ctlinedata2curtainimage)
	add_executable(${target} src/ctlinedata2curtainimage.cpp)
	target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp submodules/cpp-utils/src/gdal_utils.cpp submodules/cpp-utils/src/RamerDouglasPeucker.cpp)
	target_compile_definitions(${target} PRIVATE -DHAVE_GDAL)	
	target_include_directories(${target} PRIVATE ${GDAL_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${GDAL_LIBRARIES})
	target_link_libraries(${target} PRIVATE ticpp)
	install(TARGETS ${target} DESTINATION bin)	
else()
	message(WARNING "${target} requires GDAL - will not be built")
endif()

# Add ctlinedata2sgrid executable
set(target ctlinedata2sgrid)
add_executable(${target} src/ctlinedata2sgrid.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp) 
target_link_libraries(${target} PRIVATE ticpp)
install(TARGETS ${target} DESTINATION bin)

# Add ctlinedata2slicegrids executable
if(${WITH_GDAL} AND GDAL_FOUND) 
	set(target ctlinedata2slicegrids)
	add_executable(${target} src/ctlinedata2slicegrids.cpp)
	target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
	target_compile_definitions(${target} PRIVATE -DHAVE_GDAL)	
	target_include_directories(${target} PRIVATE ${GDAL_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${GDAL_LIBRARIES})
	install(TARGETS ${target} DESTINATION bin)
else()
	message(WARNING "${target} requires GDAL - will not be built")
endif()

# Add removelog10conductivityfromsgrid executable
set(target removelog10conductivityfromsgrid)
add_executable(${target} src/removelog10conductivityfromsgrid.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
install(TARGETS ${target} DESTINATION bin)
