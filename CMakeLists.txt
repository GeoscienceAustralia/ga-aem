cmake_minimum_required(VERSION 3.12)

# Set project name
project(ga-aem)

message(Vuild Type = ${CMAKE_BUILD_TYPE})

# Set options
option(USE_MPI "Use MPI" ON)
option(USE_NETCDF "Use NetCDF" ON)
option(USE_GDAL "Use GDAL" ON)
option(USE_PETSC "Use PETSC" ON)

message(USE_MPI = ${USE_MPI})
message(USE_NETCDF = ${USE_NETCDF})
message(USE_GDAL = ${USE_GDAL})
message(USE_PETSC = ${USE_PETSC})

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Configure find package
find_package(PkgConfig REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security -Wno-date-time")
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-unused-result -Wno-format-security -Wno-date-time")
#On GCC, even with -Wno-date-time, still get warings of the form: warning: macro "__DATE__" might prevent reproducible builds [-Wdate-time]

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	message("OpenMP was found")
endif()

# Configure FFTW
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries     (PkgConfig::FFTW)
if(FFTW_FOUND)
	message("FFTW was found")
endif()

# Configure MPI if opted for
if(${USE_MPI})
	pkg_search_module(MPI REQUIRED mpi IMPORTED_TARGET)
	if(MPI_FOUND)
		message("MPI was found")
		include_directories(PkgConfig::MPI)
		link_libraries     (PkgConfig::MPI)
	endif()
endif()

# Configure NETCDFCXX if opted for
if(${USE_NETCDF})
	pkg_search_module(NETCDFCXX REQUIRED netcdf-cxx4 IMPORTED_TARGET) 
	if(NETCDFCXX_FOUND)
		message("NETCDFCXX was found")
		include_directories(PkgConfig::NETCDFCXX)
		link_libraries     (PkgConfig::NETCDFCXX)
	    add_definitions(-DHAVE_NETCDF)
        add_definitions(-D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
	endif()
endif()

# Configure GDAL
if(${USE_GDAL})
	pkg_search_module(GDAL gdal IMPORTED_TARGET) 
	if(GDAL_FOUND)
		message("GDAL was found")
		include_directories(PkgConfig::GDAL)
		link_libraries     (PkgConfig::GDAL)
		add_definitions(-DHAVE_GDAL)
	endif()
endif()

# Configure PETSc
if(${USE_PETSC})
	pkg_search_module(PETSC PETSc REQUIRED IMPORTED_TARGET) 
	if(PETSC_FOUND)
		message("PETSC was found")
		include_directories(PkgConfig::PETSC)
		link_libraries     (PkgConfig::PETSC)
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
if(${USE_PETSC})
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

if(OpenMP_CXX_FOUND) 
	target_link_libraries(galeisbstdem.exe  PUBLIC OpenMP::OpenMP_CXX)
	target_link_libraries(galeisbsfdem.exe  PUBLIC OpenMP::OpenMP_CXX)
endif()

if(${USE_PETSC})
	if(${USE_MPI}) 
		add_executable(galeiallatonce.exe src/galeiallatonce.cpp)
		target_compile_definitions(galeiallatonce.exe PRIVATE _MPI_ENABLED OMPI_SKIP_MPICXX)
		install(TARGETS galeiallatonce.exe DESTINATION bin)
	else()
		message("galeiallatonce.exe requires PETSC and MPI")
	endif()
else()
	message("galeiallatonce.exe requires PETSC and MPI")
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
