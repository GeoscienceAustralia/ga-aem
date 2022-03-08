cmake_minimum_required(VERSION 3.12)

# Set project name
project(ga-aem)

message(STATUS Build Type = ${CMAKE_BUILD_TYPE})

# Set options
option(USE_MPI "Use MPI" ON)
option(USE_NETCDF "Use NetCDF" ON)
option(USE_GDAL "Use GDAL" ON)
option(USE_PETSC "Use PETSC" ON)

message(STATUS USE_MPI = ${USE_MPI})
message(STATUS USE_NETCDF = ${USE_NETCDF})
message(STATUS USE_GDAL = ${USE_GDAL})
message(STATUS USE_PETSC = ${USE_PETSC})

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Configure find package
find_package(PkgConfig REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")
if(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -Wno-date-time")
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-unused-result -Wno-date-time")
	#On GCC, even with -Wno-date-time, still get warings of the form: warning: macro "__DATE__" might prevent reproducible builds [-Wdate-time]
endif()

# Adding filesystem library
link_libraries(-lstdc++fs)

# Configure FFTW
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries     (PkgConfig::FFTW)
if(FFTW_FOUND)
	message(STATUS "FFTW ${FFTW_VERSION} was found")
endif()


# Configure MPI if opted for
if(${USE_MPI})
	find_package(MPI)
	if(MPI_FOUND)
		message(STATUS "MPI was found")
		include_directories(SYSTEM ${MPI_INCLUDE_PATH})
		add_definitions(-D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
		link_libraries(${MPI_C_LIBRARIES})
	endif()
endif()


# Configure GDAL
if(${USE_GDAL})
	pkg_search_module(GDAL gdal IMPORTED_TARGET) 
	if(GDAL_FOUND)
		message(STATUS "GDAL ${GDAL_VERSION} was found")
		include_directories(PkgConfig::GDAL)
		link_libraries     (PkgConfig::GDAL)
		add_definitions(-DHAVE_GDAL)
	endif()
endif()

# Configure NETCDFCXX if opted for
if(${USE_NETCDF})
	pkg_search_module(NETCDFCXX REQUIRED netcdf-cxx4 IMPORTED_TARGET) 
	if(NETCDFCXX_FOUND)
		message(STATUS "NETCDFCXX ${NETCDFCXX_VERSION} was found")
		include_directories(PkgConfig::NETCDFCXX)
		link_libraries     (PkgConfig::NETCDFCXX)
	    add_definitions(-DHAVE_NETCDF)
	endif()
endif()

# Configure PETSc (only used for galeiallatonce.exe)
if(${USE_PETSC})	
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
		message(STATUS "PETSc was NOT found -- will not build galeiallatonce.exe")
	endif()
endif()



# Add the include directories
include_directories(src/)
include_directories(submodules/cpp-utils/src/)
include_directories(submodules/csv-parser/single_include/)
include_directories(submodules/eigen/)
if(${USE_NETCDF})
	include_directories(submodules/geophysics-netcdf/src/)
	include_directories(submodules/geophysics-netcdf/submodules/marray/include/andres/)
endif()


### Add the targets

# Add the TICPP library submodule
add_subdirectory(submodules/ticpp)

# Add the cpp-utils library submodule
add_subdirectory(submodules/cpp-utils)

# Add the gatdaem1d shared library
add_library(gatdaem1d-shared SHARED src/gatdaem1d.cpp)
set_target_properties(gatdaem1d-shared PROPERTIES OUTPUT_NAME gatdaem1d)
target_link_libraries(gatdaem1d-shared PRIVATE cpp-utils-shared)
install(TARGETS gatdaem1d-shared DESTINATION bin)
install(TARGETS gatdaem1d-shared DESTINATION python)
install(TARGETS gatdaem1d-shared DESTINATION matlab)
#Install the include files
install(FILES   src/gatdaem1d.h TYPE INCLUDE)

# Add the gatdaem1d static library
add_library(gatdaem1d-static STATIC src/gatdaem1d.cpp)
set_target_properties(gatdaem1d-static PROPERTIES OUTPUT_NAME gatdaem1d)
target_link_libraries(gatdaem1d-static PRIVATE cpp-utils-static)
install(TARGETS gatdaem1d-static DESTINATION lib)

# Add ctlinedata2sgrid.exe
add_executable(ctlinedata2sgrid.exe src/ctlinedata2sgrid.cpp)
target_link_libraries(ctlinedata2sgrid.exe PRIVATE cpp-utils-static ticpp)
install(TARGETS ctlinedata2sgrid.exe DESTINATION bin)

# Add ctlinedata2slicegrids.exe
if(${USE_GDAL})
	add_executable(ctlinedata2slicegrids.exe src/ctlinedata2slicegrids.cpp)
	target_link_libraries(ctlinedata2slicegrids.exe PRIVATE cpp-utils-static)
	install(TARGETS ctlinedata2slicegrids.exe DESTINATION bin)
endif()

# Add gaforwardmodeltdem.exe
add_executable(gaforwardmodeltdem.exe src/gaforwardmodeltdem.cpp)
target_link_libraries(gaforwardmodeltdem.exe PRIVATE cpp-utils-static)
install(TARGETS gaforwardmodeltdem.exe DESTINATION bin)

# Add example_forward_model.exe
add_executable(example_forward_model.exe src/example_forward_model.cpp)
target_link_libraries(example_forward_model.exe PRIVATE cpp-utils-static)
install(TARGETS example_forward_model.exe DESTINATION bin)

# Add removelog10conductivityfromsgrid.exe
add_executable(removelog10conductivityfromsgrid.exe src/removelog10conductivityfromsgrid.cpp)
target_link_libraries(removelog10conductivityfromsgrid.exe PRIVATE cpp-utils-static)
install(TARGETS removelog10conductivityfromsgrid.exe DESTINATION bin)

# Add the pure C example executable
add_executable(example_forward_model_c.exe src/example_forward_model_c.c)
target_link_libraries(example_forward_model_c.exe PRIVATE gatdaem1d-static)
install(TARGETS example_forward_model_c.exe DESTINATION bin)

# Add the MPI executables
add_executable(galeisbstdem.exe src/galeisbstdem.cpp)
target_link_libraries(galeisbstdem.exe PRIVATE cpp-utils-static)
install(TARGETS galeisbstdem.exe DESTINATION bin)

add_executable(galeisbsfdem.exe src/galeisbsfdem.cpp)
target_link_libraries(galeisbsfdem.exe PRIVATE cpp-utils-static)
install(TARGETS galeisbsfdem.exe DESTINATION bin)

if(${USE_NETCDF} AND NETCDFCXX_FOUND)
	add_executable(garjmcmctdem.exe src/garjmcmctdem.cpp)
	target_link_libraries(garjmcmctdem.exe PRIVATE cpp-utils-static)
	install(TARGETS garjmcmctdem.exe DESTINATION bin)
else()
	message(WARNING "garjmcmctdem.exe requires NETCDF - will not be built")
endif()

# Add OpenMP 
find_package(OpenMP)
if(OpenMP_CXX_FOUND) 
	message(STATUS "OpenMP was found")
	target_link_libraries(galeisbstdem.exe PRIVATE OpenMP::OpenMP_CXX)
	target_link_libraries(galeisbsfdem.exe PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MPI_FOUND AND PETSC_FOUND) 
	add_executable(galeiallatonce.exe src/galeiallatonce.cpp)
	target_link_libraries(galeiallatonce.exe PRIVATE cpp-utils-static)
	target_include_directories(galeiallatonce.exe PRIVATE ${PETSC_INCLUDE_DIRS})
	target_link_libraries(galeiallatonce.exe PRIVATE ${PETSC_LDFLAGS})
	install(TARGETS galeiallatonce.exe DESTINATION bin)	
else()
	message(WARNING "galeiallatonce.exe requires PETSc and MPI -- will not be built")
endif()

