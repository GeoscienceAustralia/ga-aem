# This file was downloaded from https://github.com/jedbrown/cmake-modules/blob/master/FindNETCDF.cmake
# Author: Jed Brown
# Modified for ga-aem purposes

# - Find NETCDF
# Find the native NETCDF includes and library
#
#  NETCDF_INCLUDES    - where to find netcdf.h, etc
#  NETCDF_LIBRARIES   - Link these libraries when using NETCDF
#  NETCDF_FOUND       - True if NETCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C    - Just the C interface
#  NETCDF_LIBRARIES_CXX  - C++ interface, if available
#  NETCDF_LIBRARIES_F77  - Fortran 77 interface, if available
#  NETCDF_LIBRARIES_F90  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NETCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

find_path (NETCDF_INCLUDES netcdf.h
  HINTS NETCDF_DIR ENV NETCDF_DIR)


find_library (NETCDF_LIBRARIES_C       NAMES netcdf)
mark_as_advanced(NETCDF_LIBRARIES_C)

set (NETCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (NETCDF_libs "${NETCDF_LIBRARIES_C}")

get_filename_component (NETCDF_lib_dirs "${NETCDF_LIBRARIES_C}" PATH)

macro (NETCDF_check_interface lang header libs)
  if (NETCDF_${lang})
    find_path (NETCDF_INCLUDES_${lang} NAMES ${header}
      HINTS "${NETCDF_INCLUDES}" NO_DEFAULT_PATH)
    find_library (NETCDF_LIBRARIES_${lang} NAMES ${libs}
      HINTS "${NETCDF_lib_dirs}" NO_DEFAULT_PATH)
    mark_as_advanced (NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})
    if (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
      list (INSERT NETCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
    else (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
      set (NETCDF_has_interfaces "NO")
      message (STATUS "Failed to find NETCDF interface for ${lang}")
    endif (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
  endif (NETCDF_${lang})
endmacro (NETCDF_check_interface)

NETCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
NETCDF_check_interface (F77 netcdf.inc  netcdff)
NETCDF_check_interface (F90 netcdf.mod  netcdff)

set (NETCDF_LIBRARIES "${NETCDF_libs}" CACHE STRING "All NETCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NETCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES NETCDF_has_interfaces)

mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
