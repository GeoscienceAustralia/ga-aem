if(MSVC)
	message(STATUS "An MSVC Import Library will be created for FFTW ...")
	#For Visual studio/MSVC make an import library from the downloaded FFTW dll (something to do with it being built using MinGW)
	#and must add MSVC_FFTW_IMPORTLIB as a dependency (see add_dependencies()) for all targets that require FFTW in MSVC
	set(FFTW_LINK_LIBRARIES ${FFTW_DIR}/libfftw3-3.lib)
	add_custom_command(OUTPUT ${FFTW_DIR}/libfftw3-3.lib
		COMMAND CMD /C lib.exe /machine:x64 /def:libfftw3-3.def /out:libfftw3-3.lib
		WORKING_DIRECTORY ${FFTW_DIR}
		COMMENT "Creating MSVC FFTW import library"
		VERBATIM)
	add_custom_target(MSVC_FFTW_IMPORTLIB DEPENDS ${FFTW_DIR}/libfftw3-3.lib)
endif()
