#For Visual studio/MSVC make an import library from the downloaded FFTW dll (something to do with it being built using MinGW)
#and must add MSVC_FFTW_IMPORTLIB as a dependency (see add_dependencies()) for all targets that require FFTW in MSVC
set(fftwimportlib libfftw3-3)
set(fftwimportlib_wdpath ${CMAKE_CURRENT_BINARY_DIR})
set(fftwimportlib_dllpath ${FFTW_DIR}/${fftwimportlib}.dll)
set(fftwimportlib_defpath ${FFTW_DIR}/${fftwimportlib}.def)
set(fftwimportlib_libpath ${fftwimportlib_wdpath}/${fftwimportlib}.lib)

set(FFTW_LINK_LIBRARIES ${fftwimportlib_libpath})

message(STATUS "A MSVC FFTW import library ${fftwimportlib_libpath} will be generated and used")
add_custom_command(OUTPUT ${fftwimportlib_libpath}
	COMMAND CMD /C lib.exe /machine:x64 /def:${fftwimportlib_defpath} /out:${fftwimportlib_libpath}
	WORKING_DIRECTORY ${fftwimportlib_wdpath}
	COMMENT "Creating MSVC import library ${fftwimportlib_libpath}"
	VERBATIM)

add_custom_target(generate_msvc_import_library DEPENDS ${fftwimportlib_libpath})

set(fftwimportlib)
set(fftwimportlib_wdpath)
set(fftwimportlib_dllpath)
set(fftwimportlib_defpath)
set(fftwimportlib_libpath)

