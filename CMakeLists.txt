cmake_minimum_required(VERSION 3.12)

cmake_policy(SET CMP0076 NEW)

function(reportvar var)
  message(STATUS "${ARGVN}${var}=${${var}}")
endfunction()

function(reportallvars)
	get_cmake_property(_variableNames VARIABLES)
	list (SORT _variableNames)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endfunction()

# Set project name
project(ga-aem)

reportvar(CMAKE_SYSTEM_NAME)
reportvar(CMAKE_BUILD_TYPE)

# Set options
option(WITH_MPI "WITH_MPI" ON)
option(WITH_NETCDF "WITH_NETCDF" ON)
option(WITH_GDAL "WITH_GDAL" ON)
option(WITH_PETSC "WITH_PETSC" ON)

reportvar(WITH_MPI)
reportvar(WITH_NETCDF)
reportvar(WITH_GDAL)
reportvar(WITH_PETSC)

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

# Find PkgConfig
find_package(PkgConfig QUIET)

# Add the project include directories
include_directories(src/)
include_directories(submodules/cpp-utils/src/)
include_directories(submodules/csv-parser/single_include/)
include_directories(submodules/eigen/)

# Configure FFTW
if(PkgConfig_FOUND)
	message(STATUS "FFTW package was found")
	pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
endif()
if(NOT FFTW_FOUND)
	message(STATUS "FFTW package was not found - manual setup")
	set(FFTW_DIR $ENV{FFTW_DIR})
	find_path(FFTW_INCLUDE_DIRS NAMES "fftw3.h" PATHS ${FFTW_DIR} PATH_SUFFIXES "include" NO_DEFAULT_PATH)
	find_library(FFTW_LINK_LIBRARIES "fftw3-3" libfftw3-3 PATHS ${FFTW_DIR} PATH_SUFFIXES "lib" "lib64" NO_DEFAULT_PATH)
	set(FFTW_FOUND 1)
endif()

if(FFTW_FOUND)
	reportvar(FFTW_VERSION "	")
	reportvar(FFTW_INCLUDE_DIRS "	")
	reportvar(FFTW_LINK_LIBRARIES "	")
	reportvar(FFTW_LIBDIR "	")
	reportvar(FFTW_LIBRARIES "	")
	reportvar(FFTW_LDFLAGS "	")
endif()

#$find_package(FFTW IMPORTED_TARGET)
#if(FFTW_FOUND)
#	message(STATUS "FFTW ${FFTW_VERSION} was found")
#else()
#	set (FFTW_INCLUDE_DIRS ${FFTW_DIR})
#	set (FFTW_LIBRARIES ${FFTW_DIR})
#endif()


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
		reportvar(GDAL_INCLUDE_DIRS "	")
		reportvar(GDAL_LIBRARIES "	")
	endif()
endif()

# Configure PETSc (only used for galeiallatonce)
if(${WITH_PETSC})	
	SET(PETSC_FOUND 1)
	set(PETSC_DIR $ENV{PETSC_DIR})
	set(PETSC_INCLUDE_DIRS $ENV{PETSC_INCLUDE_DIRS})
	set(PETSC_LIB_DIR $ENV{PETSC_LIB_DIR})
	set(PETSC_LIBRARIES $ENV{PETSC_LIBRARIES})
	set(PETSC_LDFLAGS $ENV{PETSC_LDFLAGS})

	foreach (lib_name ${PETSC_LIBRARIES})
		set(full_lib_name ${PETSC_LIB_DIR}/${lib_name})
		cmake_path(SET full_lib_name NORMALIZE "${full_lib_name}")
		list(APPEND PETSC_LINK_LIBRARIES ${full_lib_name})
	endforeach()

	message(STATUS "PETSC was found")	
	reportvar(PETSC_DIR "	")	
	reportvar(PETSC_INCLUDE_DIRS "	")	
	reportvar(PETSC_LIB_DIR "	")	
	reportvar(PETSC_LIBRARIES "	")
	reportvar(PETSC_LINK_LIBRARIES "	")
	reportvar(PETSC_LDFLAGS "	")
endif()

if(${WITH_PETSC_ZZZ})
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
		message(STATUS "PETSc was NOT found")
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


# Configure NETCDF if opted for
if(${WITH_NETCDF})
	pkg_search_module(NETCDFCXX netcdf-cxx4 IMPORTED_TARGET) 
	if(NETCDFCXX_FOUND)
		message(STATUS "NETCDFCXX package was found")
		set(NETCDFCXX_PACKAGE_FOUND 1)
	else()
		message(STATUS "NETCDFCXX package was not found - building our own")
		set(NETCDF_DIR $ENV{NETCDF_DIR})
		set(NETCDF_INCLUDE_DIR $ENV{NETCDF_INCLUDE_DIR})
		set(NETCDF_LIB_DIR $ENV{NETCDF_LIB_DIR})
		set(NETCDF_LIBRARIES $ENV{NETCDF_LIBRARIES})

		foreach (lib_name ${NETCDF_LIBRARIES})
			set(full_lib_name ${NETCDF_LIB_DIR}/${lib_name})
			cmake_path(SET full_lib_name NORMALIZE "${full_lib_name}")
			list(APPEND NETCDF_LINK_LIBRARIES ${full_lib_name})
		endforeach()

		reportvar(NETCDF_DIR "	")
		reportvar(NETCDF_INCLUDE_DIR "	")
		reportvar(NETCDF_LIB_DIR "	")
		reportvar(NETCDF_LIBRARIES "	")
		reportvar(NETCDF_LINK_LIBRARIES "	")

		#set(NETCDF_LINK_LIBRARIES 
		#	${NETCDF_LIB_DIR}/netcdf.lib 
		#	${NETCDF_LIB_DIR}/hdf5.lib 
		#	${NETCDF_LIB_DIR}/hdf5_hl.lib 
		#	${NETCDF_LIB_DIR}/libcurl_imp.lib)
		set(NETCDFCXX_SRC_DIR submodules/netcdf-cxx4/cxx4)
		set(NETCDFCXX_INCLUDE_DIRS ${NETCDFCXX_SRC_DIR} ${NETCDF_INCLUDE_DIR})
		set(NETCDFCXX_PACKAGE_FOUND 0)


		file(GLOB NC_CXX_HEADERS ${NETCDFCXX_SRC_DIR}/*.h)
		file(GLOB NC_CXX_SOURCES ${NETCDFCXX_SRC_DIR}/nc*.cpp)
		add_library(libnetcdfcxx STATIC ${NC_CXX_SOURCES})
		target_include_directories(libnetcdfcxx PUBLIC ${NETCDFCXX_INCLUDE_DIRS})
		target_link_libraries(libnetcdfcxx PUBLIC ${NETCDF_LINK_LIBRARIES}) 
		#include_directories(libnetcdfcxx INTERFACE ${netCDF_INCLUDE_DIR})
		#include_directories(submodules/geophysics-netcdf/src/)
		#include_directories(submodules/geophysics-netcdf/submodules/marray/include/andres/)
		reportvar(NETCDFCXX_INCLUDE_DIRS "	")
		reportvar(NETCDFCXX_LIBRARY_DIRS "	")	
		reportvar(NETCDFCXX_LIBRARIES "	")
		reportvar(NETCDFCXX_LINK_LIBRARIES "	")
		reportvar(NETCDFCXX_LD_FLAGS "	")
	endif()

endif()

### Add the targets
# Add the NetCDF CXX library
if(${WITH_NETCDFZZZ})	

endif()

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
target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	
install(TARGETS ${target} DESTINATION lib)

# Add the gatdaem1d shared library
set(target gatdaem1d-shared)
add_library(${target} SHARED src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
target_sources(${target} PRIVATE submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	
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
target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	
install(TARGETS ${target} DESTINATION bin)

# Add example_forward_model executable
set(target example_forward_model)
add_executable(${target} src/example_forward_model.cpp)
target_sources(${target} PUBLIC submodules/cpp-utils/src/file_utils.cpp submodules/cpp-utils/src/general_utils.cpp)
target_include_directories(${target} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	
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
target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})
if(${WITH_MPI})
	target_compile_definitions(${target} PRIVATE -D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
	target_include_directories(${target} PRIVATE SYSTEM ${MPI_INCLUDE_PATH})
	target_link_libraries(${target} PRIVATE ${MPI_C_LIBRARIES})
endif()
if(${WITH_NETCDF})
	target_compile_definitions(${target} PRIVATE -DHAVE_NETCDF)
	target_include_directories(${target} PRIVATE submodules/geophysics-netcdf/src/)
	target_include_directories(${target} PRIVATE submodules/geophysics-netcdf/submodules/marray/include/andres/)	
	target_include_directories(${target} PRIVATE ${NETCDFCXX_INCLUDE_DIRS})

	target_link_libraries(${target} PRIVATE ${NETCDFCXX_LINK_LIBRARIES})
	if(NOT ${NETCDFCXX_PACKAGE_FOUND})
		target_link_libraries(${target} PRIVATE libnetcdfcxx)
	endif()
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
	target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	

	target_compile_definitions(${target} PRIVATE -DHAVE_NETCDF)
	target_include_directories(${target} PRIVATE ${NETCDFCXX_INCLUDE_DIRS})
	target_link_libraries(${target} PRIVATE ${NETCDFCXX_LINK_LIBRARIES})
	if(NOT ${NETCDFCXX_PACKAGE_FOUND})
		target_link_libraries(${target} PRIVATE libnetcdfcxx)
	endif()	

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
	target_link_libraries(${target} PRIVATE ${FFTW_LINK_LIBRARIES})	

	target_compile_definitions(${target} PRIVATE -D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
	target_include_directories(${target} PRIVATE SYSTEM ${MPI_INCLUDE_PATH})
	target_link_libraries(${target} PRIVATE ${MPI_C_LIBRARIES})

	target_include_directories(${target} PRIVATE ${PETSC_INCLUDE_DIRS})
	#target_link_libraries(${target} PRIVATE ${PETSC_LINK_LIBRARIES})	
	target_link_directories(${target} PRIVATE ${PETSC_LIB_DIR})
	target_link_libraries(${target} PRIVATE ${PETSC_LIBRARIES})
	target_link_options(${target} PRIVATE ${PETSC_LDFLAGS})

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
