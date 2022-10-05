message(STATUS "Running cmake post install script")

if(UNIX)
	message(STATUS "Copying shared library libgatdaem1d.so to matlab and python directories")
	configure_file("${CMAKE_INSTALL_PREFIX}/bin/libgatdaem1d.so" "${CMAKE_INSTALL_PREFIX}/python/gatdaem1d/gatdaem1d.so" COPYONLY)	
	configure_file("${CMAKE_INSTALL_PREFIX}/bin/libgatdaem1d.so" "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.so" COPYONLY)	
else()
	message(STATUS "Copying shared library libgatdaem1d.dll to matlab and python directories")
	configure_file("${CMAKE_INSTALL_PREFIX}/bin/gatdaem1d.dll" "${CMAKE_INSTALL_PREFIX}/python/gatdaem1d/gatdaem1d.dll" COPYONLY)
	configure_file("${CMAKE_INSTALL_PREFIX}/bin/gatdaem1d.dll" "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.mexw64" COPYONLY)
endif()

