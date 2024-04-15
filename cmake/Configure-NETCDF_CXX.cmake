pkg_search_module(NETCDFCXX IMPORTED_TARGET netcdf-cxx4)
if(NETCDFCXX_FOUND)
	message(STATUS "NETCDF C++ libraries were found")
	add_library(NETCDF::CXX INTERFACE IMPORTED)
	set_target_properties(NETCDF::CXX PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${NETCDFCXX_INCLUDE_DIRS}"
		INTERFACE_LINK_DIRECTORIES "${NETCDFCXX_LIBRARY_DIRS}"
		INTERFACE_LINK_LIBRARIES "${NETCDFCXX_LINK_LIBRARIES}"
		INTERFACE_LINK_OPTIONS "${NETCDFCXX_LINK_OPTIONS}"
	)
else()
	message(STATUS "NETCDF C++ libraries were NOT found -- we will build our own locally")
	set(NETCDF_DIR  $ENV{NETCDF_DIR})
	set(NETCDF_ROOT $ENV{NETCDF_DIR})
	find_package(NETCDF REQUIRED)
	if(NOT NETCDF_FOUND)
		message(FATAL_ERROR "The NETCDF C Libraries could not be found")
	endif()
	if(NOT NETCDFCXX_SRC_DIR)
		message(FATAL_ERROR "The NETCDFCXX_SRC_DIR variable was not set")
	else()
		set(NETCDFCXX_INCLUDE_DIRS ${NETCDFCXX_SRC_DIR})
		file(GLOB NC_CXX_SOURCES ${NETCDFCXX_SRC_DIR}/nc*.cpp)

		set(target netcdf-cxx4)
		add_library(${target} STATIC ${NC_CXX_SOURCES})
		set_target_properties(${target} PROPERTIES OUTPUT_NAME ${target})
		target_include_directories(${target} PUBLIC ${NETCDF_INCLUDES})
		target_include_directories(${target} PUBLIC ${NETCDFCXX_INCLUDE_DIRS})
		target_link_libraries(${target} PUBLIC ${NETCDF_LIBRARIES_C})
		add_library(NETCDF::CXX ALIAS netcdf-cxx4)

		set(NETCDFCXX_FOUND TRUE)
	endif()
endif()

if(NETCDFCXX_FOUND)
	#reportvar(NETCDFCXX_INCLUDE_DIRS)
	#reportvar(NETCDFCXX_LIBRARY_DIRS)
	#reportvar(NETCDFCXX_LINK_LIBRARIES)
	#reportvar(NETCDFCXX_LINK_OPTIONS)
endif()
