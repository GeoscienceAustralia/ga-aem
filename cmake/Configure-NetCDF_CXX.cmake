message(STATUS "\nChecking for the NETCDF C++ libraries")
pkg_search_module(NETCDFCXX IMPORTED_TARGET netcdf-cxx4)
if(NETCDFCXX_FOUND)
	message(STATUS "NETCDF C++ libraries were found")
	message(STATUS "Creating NETCDF::CXX INTERFACE")
	add_library(NETCDF::CXX INTERFACE IMPORTED)
	set_target_properties(NETCDF::CXX PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${NETCDFCXX_INCLUDE_DIRS}"
		INTERFACE_LINK_DIRECTORIES "${NETCDFCXX_LIBRARY_DIRS}"
		INTERFACE_LINK_LIBRARIES "${NETCDFCXX_LINK_LIBRARIES}"
		INTERFACE_LINK_OPTIONS "${NETCDFCXX_LINK_OPTIONS}"
	)
else()
	message(STATUS "NETCDF C++ libraries were NOT found -- we will build our own locally")
	include(cmake/Configure-NetCDF_C.cmake)
	set(NETCDFCXX_SRC_DIR submodules/netcdf-cxx4/cxx4)
	set(NETCDFCXX_INCLUDE_DIRS ${NETCDFCXX_SRC_DIR})
	file(GLOB NC_CXX_SOURCES ${NETCDFCXX_SRC_DIR}/nc*.cpp)

	set(target netcdf_c++4)
	add_library(${target} STATIC ${NC_CXX_SOURCES})
	add_dependencies(${target} NETCDF::C)
	set_target_properties(${target} PROPERTIES OUTPUT_NAME ${target})
	target_include_directories(${target} PUBLIC ${NETCDFCXX_INCLUDE_DIRS})
	target_link_libraries(${target} PUBLIC NETCDF::C)
	set(NETCDFCXX_FOUND TRUE)

	add_library(NETCDF::CXX ALIAS netcdf_c++4)
endif()

if(NETCDFCXX_FOUND)
	reportvar(NETCDFCXX_INCLUDE_DIRS)
	reportvar(NETCDFCXX_LIBRARY_DIRS)
	reportvar(NETCDFCXX_LINK_LIBRARIES)
	reportvar(NETCDFCXX_LINK_OPTIONS)
endif()
