## Do NOT edit this file unless you really know what you are doing
cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0076 NEW)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_COLOR_MAKEFILE ON)

## Set project name
project(ga-aem VERSION 2.0 DESCRIPTION "Airborne electromagnetics forward modelling and inversion" LANGUAGES CXX C)
message(STATUS "Configuring ga-aem")

include(cmake/Message-Functions.cmake)
include(cmake/Check-Functions.cmake)

# Set options
option(WITH_FFTW "Build with FFTW support" ON)
option(WITH_MPI  "Build with MPI support" ON)
option(WITH_NETCDF "Build with NetCDF support" ON)
option(WITH_GDAL "Build with GDAL support" ON)
option(WITH_PETSC "Build with PETSC support" ON)

reportvar(WITH_FFTW)
reportvar(WITH_MPI)
reportvar(WITH_NETCDF)
reportvar(WITH_GDAL)
reportvar(WITH_PETSC)

reportvar(CMAKE_VERSION)
reportvar(CMAKE_SYSTEM_NAME)
reportvar(CMAKE_BUILD_TYPE)

set(CMAKE_EXECUTABLE_SUFFIX_C ".exe")
set(CMAKE_EXECUTABLE_SUFFIX_CXX ".exe")

# Include the project level options and flags
include(cmake/Typical-Compiler-Options.cmake)

# Find Configure the external packages
#This is where the NetCDF C++ src code is located if it is not installed and has to be built from source
set(NETCDFCXX_SRC_DIR submodules/netcdf-cxx4/cxx4)
include(cmake/Configure-External-Packages.cmake)

#######################################
### Add the build targets

# Add the project include directories, which don't have their own CMakelists.txt
include_directories(src/)
include_directories(submodules/eigen)

# Add the cpp-utils submodule
add_subdirectory(submodules/cpp-utils EXCLUDE_FROM_ALL)

# Add the geophysics-netcdf submodule
if(${WITH_NETCDF})
	add_subdirectory(submodules/geophysics-netcdf EXCLUDE_FROM_ALL)
endif()

# Add the TICPP library submodule
add_subdirectory(submodules/ticpp EXCLUDE_FROM_ALL)

# Add the gatdaem1d static library
set(target gatdaem1d-static)
add_library(${target} STATIC EXCLUDE_FROM_ALL src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d-static)
set_target_properties(${target} PROPERTIES PREFIX "")
set_target_properties(${target} PROPERTIES PUBLIC_HEADER src/gatdaem1d.h)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
install(TARGETS ${target} OPTIONAL PUBLIC_HEADER DESTINATION include OPTIONAL)

# Platform differences for Matlab and Python targets
if(WIN32)
	set(matlab_suffix ".mexw64")
	set(artefact_type "RUNTIME")
else()
	set(matlab_suffix ".so")
	set(artefact_type "LIBRARY")
endif()

# Add the gatdaem1d shared library
set(target gatdaem1d-shared)
add_library(${target} SHARED EXCLUDE_FROM_ALL src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
set_target_properties(${target} PROPERTIES PREFIX "")
set_target_properties(${target} PROPERTIES PUBLIC_HEADER src/gatdaem1d.h)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
install(TARGETS ${target} ${artefact_type} OPTIONAL PUBLIC_HEADER DESTINATION include OPTIONAL)

# Python bindings target - same as gatdaem1d-shared but different location
set(target python-bindings)
add_library(${target} SHARED src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
set_target_properties(${target} PROPERTIES PREFIX "")
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
install(TARGETS ${target} ${artefact_type} DESTINATION python/gatdaem1d OPTIONAL)

# Matlab bindings target - same as gatdaem1d-shared but different suffix and location
set(target matlab-bindings)
add_library(${target} SHARED src/gatdaem1d.cpp)
set_target_properties(${target} PROPERTIES OUTPUT_NAME gatdaem1d)
set_target_properties(${target} PROPERTIES PREFIX "")
set_target_properties(${target} PROPERTIES SUFFIX ${matlab_suffix})
set_target_properties(${target} PROPERTIES PUBLIC_HEADER src/gatdaem1d.h)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
install(TARGETS ${target} ${artefact_type} DESTINATION matlab/bin OPTIONAL PUBLIC_HEADER DESTINATION matlab/gatdaem1d_functions OPTIONAL)

# Add gaforwardmodeltdem executable
set(target gaforwardmodeltdem)
add_executable(${target} src/gaforwardmodeltdem.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add example_forward_model executable
set(target example_forward_model)
add_executable(${target} src/example_forward_model.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
# There is no need to install this example as its just for developers
#install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add the pure C example executable
set(target example_forward_model_c)
add_executable(${target} src/example_forward_model_c.c)
target_link_libraries(${target} PRIVATE gatdaem1d-static)
target_link_libraries(${target} PRIVATE cpp-utils)
# There is no need to install this example as its just for developers
#install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add galeisbstdem executable
set(target galeisbstdem)
add_executable(${target} src/galeisbstdem.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
if(${WITH_MPI})
	target_compile_definitions(${target} PRIVATE ENABLE_MPI OMPI_SKIP_MPICXX)
endif()
if(${WITH_NETCDF})
	target_compile_definitions(${target} PRIVATE HAVE_NETCDF)
	target_link_libraries(${target} PUBLIC geophysics-netcdf)
endif()
if(OpenMP_CXX_FOUND)
	target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
endif()
install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add galeisbstdem-nompi executable
# This is simply an executable without any MPI functionality or dependency designed for use by people who are not allowed to have MPI installed on their systems
set(target galeisbstdem-nompi)
add_executable(${target} src/galeisbstdem.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE FFTW::FFTW)
if(${WITH_NETCDF})
	target_compile_definitions(${target} PRIVATE HAVE_NETCDF)
	target_link_libraries(${target} PUBLIC geophysics-netcdf)
endif()
if(OpenMP_CXX_FOUND)
	target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
endif()
install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add garjmcmctdem executable
set(target garjmcmctdem)
if(${WITH_NETCDF})
	add_executable(${target} src/garjmcmctdem.cpp)
	target_link_libraries(${target} PRIVATE cpp-utils)
	target_link_libraries(${target} PRIVATE FFTW::FFTW)
	target_compile_definitions(${target} PRIVATE HAVE_NETCDF)
	target_link_libraries(${target} PUBLIC geophysics-netcdf)
	if(${WITH_MPI})
		target_compile_definitions(${target} PRIVATE ENABLE_MPI OMPI_SKIP_MPICXX)
	endif()
	install(TARGETS ${target} DESTINATION bin OPTIONAL)
else()
	message(WARNING "${target} requires NETCDF - will not be built")
endif()

# Add galeiallatonce executable
set(target galeiallatonce)
if(MPI_C_FOUND AND PETSC_FOUND)
	add_executable(${target} src/galeiallatonce.cpp)
	target_link_libraries(${target} PRIVATE cpp-utils)
	target_link_libraries(${target} PRIVATE FFTW::FFTW)
	target_compile_definitions(${target} PRIVATE ENABLE_MPI OMPI_SKIP_MPICXX)
	target_link_libraries(${target} PRIVATE PETSC::PETSC)
	install(TARGETS ${target} DESTINATION bin OPTIONAL)
else()
	message(WARNING "${target} requires PETSc and MPI -- will not be built")
endif()

# Add galeisbsfdem executable
### This frequency-domain inversion code is currently out of action
### set(target galeisbsfdem)
### add_executable(${target} src/galeisbsfdem.cpp)
### target_link_libraries(${target} PRIVATE cpp-utils)
### if(${WITH_MPI})
### 		target_compile_definitions(${target} PRIVATE ENABLE_MPI OMPI_SKIP_MPICXX)
### endif()
### if(OpenMP_CXX_FOUND)
### 	target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
### endif()
### install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add ctlinedata2sgrid executable
set(target ctlinedata2sgrid)
add_executable(${target} src/ctlinedata2sgrid.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
target_link_libraries(${target} PRIVATE ticpp::ticpp)
install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Add ctlinedata2slicegrids executable
set(target ctlinedata2slicegrids)
if(${WITH_GDAL} AND GDAL_FOUND)
	add_executable(${target} src/ctlinedata2slicegrids.cpp)
	target_link_libraries(${target} PRIVATE cpp-utils)
	install(TARGETS ${target} DESTINATION bin OPTIONAL)
else()
	message(WARNING "${target} requires GDAL - will not be built")
endif()

# Add ctlinedata2georefimage executable
set(target ctlinedata2georefimage)
if (MSVC)
	# Only on Windows with MSVC
	add_executable(${target} src/ctlinedata2georefimage.cpp)
	target_link_libraries(${target} PRIVATE cpp-utils)
	install(TARGETS ${target} DESTINATION bin OPTIONAL)
endif()

# Add ctlinedata2curtainimage executable
set(target ctlinedata2curtainimage)
if (MSVC)
	# Only on Windows with MSVC
	if(${WITH_GDAL} AND GDAL_FOUND)
		add_executable(${target} src/ctlinedata2curtainimage.cpp)
		target_link_libraries(${target} PRIVATE cpp-utils)
		target_link_libraries(${target} PRIVATE ticpp::ticpp)
		install(TARGETS ${target} DESTINATION bin OPTIONAL)
	else()
		message(WARNING "${target} requires GDAL - will not be built")
	endif()
endif()

# Add removelog10conductivityfromsgrid executable
set(target removelog10conductivityfromsgrid)
add_executable(${target} src/removelog10conductivityfromsgrid.cpp)
target_link_libraries(${target} PRIVATE cpp-utils)
install(TARGETS ${target} DESTINATION bin OPTIONAL)

# Install the documentation
install(DIRECTORY docs/ DESTINATION docs MESSAGE_NEVER)

# Install the examples
install(DIRECTORY examples/ DESTINATION examples
		MESSAGE_NEVER
		PATTERN "*.sh" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)

# Install the python library and examples and code, excluding rubbish from the development environment
install(DIRECTORY python/ DESTINATION python
		MESSAGE_NEVER
		PATTERN "*.ilk" EXCLUDE
		PATTERN "*.pdb" EXCLUDE
		PATTERN "*.ipdb" EXCLUDE
		PATTERN "*.iobj" EXCLUDE
		PATTERN "*.dll" EXCLUDE
		PATTERN "__pycache__" EXCLUDE
		PATTERN "*.sh" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)

# Install the matlab library and examples and functions, excluding rubbish from the development environment
install(DIRECTORY matlab/ DESTINATION matlab
		MESSAGE_NEVER
		PATTERN "bin" EXCLUDE)

# Install any Unix shell scripts or Windows batch files
if(UNIX)
	install(DIRECTORY scripts/ DESTINATION scripts
			MESSAGE_NEVER
			FILES_MATCHING PATTERN "*.sh"
			PATTERN "*.sh" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)
else()
	install(DIRECTORY scripts/ DESTINATION scripts
			MESSAGE_NEVER
			FILES_MATCHING PATTERN "*.bat"
			PATTERN "*.bat" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)
endif()

# Install certain root directory files
install(FILES README.md DESTINATION .)
install(FILES README-Dependencies.md DESTINATION .)
install(FILES COPYRIGHT.txt DESTINATION .)
install(FILES LICENCE.txt DESTINATION .)

