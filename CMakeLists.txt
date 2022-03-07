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

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	message(STATUS "OpenMP was found")
endif()

# Configure MPI if opted for
if(${USE_MPI})
	find_package(MPI)
	if(MPI_FOUND)
		message(STATUS "MPI was found")
		include_directories(SYSTEM ${MPI_INCLUDE_PATH})
		add_definitions(-DOMPI_SKIP_MPICXX)
	endif()

	#pkg_search_module(MPI REQUIRED mpi IMPORTED_TARGET)
	#if(MPI_FOUND)
		#message("MPI was found")
		#include_directories(PkgConfig::MPI)
		#link_libraries     (PkgConfig::MPI)
	#endif()
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

# Configure PETSc
if(${USE_PETSC})	
	pkg_search_module(PETSC PETSc IMPORTED_TARGET) 
	if(PETSC_FOUND)		
		message(STATUS "PETSc ${PETSc_VERSION} was found")		
		#message("${PETSC_INCLUDE_DIRS}")
		#message("${PETSC_LIBDIR}")
		#message("${PETSC_LDFLAGS}")
	else()
		# If PETSc is not set up properly try a manual approach
		find_package(PETSc)	
		if(PETSC_FOUND)		
			set(PETSC_LIBDIR /apps/petsc/3.12.2/lib/ompi3/Intel)
			set(PETSC_LDFLAGS "${PETSC_LDFLAGS} -L/apps/intel-ct/2019.5.281/mkl/lib/intel64 -L/apps/intel-ct/2019.5.281/mkl/../compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/8 -L/apps/netcdf/4.7.1p/lib -L/apps/hdf5/1.10.5p/lib -L/apps/fftw3/3.3.8/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lsuperlu -lsuperlu_dist -lml -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lirc -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lnetcdf -lparmetis -lmetis -ltriangle -lm -lz -lX11 -lctetgen -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lrt -lquadmath -lstdc++ -ldl")
			#message(STATUS "PETSc ${PETSC_VERSION} was found")		
			#message(STATUS "PETSc VERSION ${PETSC_VERSION}")		
			#message(STATUS "PETSC_LIBDIR=${PETSC_LIBDIR}")	
			#message(STATUS "PETSC_LDFLAGS=${PETSC_LDFLAGS}")	
		else()
			message(STATUS "PETSc was NOT found")
		endif()
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

# Add the TICPP library submodule
add_subdirectory(submodules/ticpp)

# Add the cpp-utils library submodule
add_subdirectory(submodules/cpp-utils)

# Add the gatdaem1d shared library
add_library(gatdaem1d-shared SHARED src/gatdaem1d.cpp)
target_link_libraries(gatdaem1d-shared PUBLIC cpp-utils-shared)
set_target_properties(gatdaem1d-shared PROPERTIES OUTPUT_NAME gatdaem1d)

# Add the gatdaem1d static library
add_library(gatdaem1d-static STATIC src/gatdaem1d.cpp)
target_link_libraries(gatdaem1d-static PUBLIC cpp-utils-static)
set_target_properties(gatdaem1d-static PROPERTIES OUTPUT_NAME gatdaem1d)

# Link cpp-utils against all targets
link_libraries(cpp-utils-static)

# Add the NON-MPI executables
add_executable(ctlinedata2sgrid.exe src/ctlinedata2sgrid.cpp)
target_link_libraries(ctlinedata2sgrid.exe PUBLIC ticpp)
if(${USE_GDAL})
	add_executable(ctlinedata2slicegrids.exe src/ctlinedata2slicegrids.cpp)
	install(TARGETS ctlinedata2slicegrids.exe DESTINATION bin)
endif()

add_executable(gaforwardmodeltdem.exe src/gaforwardmodeltdem.cpp)
add_executable(example_forward_model.exe src/example_forward_model.cpp)
add_executable(removelog10conductivityfromsgrid.exe src/removelog10conductivityfromsgrid.cpp)

# Add the NON-MPI pure C example executable
add_executable(example_forward_model_c.exe src/example_forward_model_c.c)
target_link_libraries(example_forward_model_c.exe PUBLIC gatdaem1d-static)

# Add the MPI executables
add_executable(galeisbstdem.exe src/galeisbstdem.cpp)
if(${USE_NETCDF})
	add_executable(garjmcmctdem.exe src/garjmcmctdem.cpp)
endif()
add_executable(galeisbsfdem.exe src/galeisbsfdem.cpp)


# Add OpenMP 
if(OpenMP_CXX_FOUND) 
	target_link_libraries(galeisbstdem.exe  PUBLIC OpenMP::OpenMP_CXX)
	target_link_libraries(galeisbsfdem.exe  PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MPI_FOUND AND PETSC_FOUND) 
	add_executable(galeiallatonce.exe src/galeiallatonce.cpp)
	target_link_libraries(galeiallatonce.exe PUBLIC  ${MPI_C_LIBRARIES})

	target_include_directories(galeiallatonce.exe PRIVATE ${PETSC_INCLUDE_DIRS})
	#target_link_libraries(galeiallatonce.exe PRIVATE "-L${PETSC_LIBDIR} ${PETSC_LDFLAGS}")
	target_link_libraries(galeiallatonce.exe PRIVATE ${PETSC_LDFLAGS})
	install(TARGETS galeiallatonce.exe DESTINATION bin)	
else()
	message(WARNING "galeiallatonce.exe requires PETSc and MPI -- will not be built")
endif()

if(${USE_MPI} AND MPI_FOUND) 
	target_link_libraries(galeisbstdem.exe PUBLIC ${MPI_C_LIBRARIES})
	target_link_libraries(galeisbsfdem.exe PUBLIC ${MPI_C_LIBRARIES})
	target_link_libraries(garjmcmctdem.exe PUBLIC ${MPI_C_LIBRARIES})	
endif()


install(TARGETS galeisbstdem.exe DESTINATION bin)

if(${USE_NETCDF})
	install(TARGETS garjmcmctdem.exe DESTINATION bin)
endif()

install(TARGETS galeisbsfdem.exe DESTINATION bin)
install(TARGETS ctlinedata2sgrid.exe DESTINATION bin)
install(TARGETS example_forward_model.exe DESTINATION bin)
install(TARGETS example_forward_model_c.exe DESTINATION bin)
install(TARGETS gaforwardmodeltdem.exe DESTINATION bin)
install(TARGETS removelog10conductivityfromsgrid.exe DESTINATION bin)

install(TARGETS gatdaem1d-shared DESTINATION bin)
install(TARGETS gatdaem1d-shared DESTINATION python)
install(TARGETS gatdaem1d-shared DESTINATION matlab)
install(TARGETS gatdaem1d-static DESTINATION lib)
install(FILES   src/gatdaem1d.h TYPE INCLUDE)
