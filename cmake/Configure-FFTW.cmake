set(FFTW_DIR $ENV{FFTW_DIR})

pkg_search_module(FFTW IMPORTED_TARGET fftw3)
if(NOT FFTW_FOUND)
	message(STATUS "FFTW was not found by pkg_search_module()")
endif()

if(NOT FFTW_FOUND)
	if (MSVC)
		message(STATUS "FFTW -- resorting to manual setup")
		find_path(FFTW_INCLUDE_DIRS NAMES "fftw3.h" PATHS ${FFTW_DIR} PATH_SUFFIXES "include" NO_DEFAULT_PATH)
		check_existance(FFTW_INCLUDE_DIRS status)
		if(NOT status)
			message(FATAL_ERROR "The FFTW include directories were not found -- make sure FFTW_DIR variable is set in your environment")
		endif()
		include(cmake/Create-MSVC-FFTW-Import-Library.cmake)
		set(FFTW_FOUND TRUE)
	endif()
endif()

if(FFTW_FOUND)
	set(FFTW_LINK_LIBRARIES ${FFTW_LINK_LIBRARIES})
	set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDEDIR})
	if (NOT TARGET FFTW::FFTW)
		add_library(FFTW::FFTW INTERFACE IMPORTED)
		set_target_properties(FFTW::FFTW PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
			INTERFACE_LINK_LIBRARIES "${FFTW_LINK_LIBRARIES}"
		)
		if(MSVC)
			add_dependencies(FFTW::FFTW generate_msvc_import_library)
		endif()
	endif ()
endif ()


