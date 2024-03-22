message(STATUS "Running cmake post install script")

if(UNIX)
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libgatdaem1d.so")
		message(STATUS "Copying shared library libgatdaem1d.so to matlab and python directories")
		configure_file("${CMAKE_INSTALL_PREFIX}/lib/libgatdaem1d.so" "${CMAKE_INSTALL_PREFIX}/python/gatdaem1d/gatdaem1d.so" COPYONLY)	
		configure_file("${CMAKE_INSTALL_PREFIX}/lib/libgatdaem1d.so" "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.so" COPYONLY) 
	endif()
else()
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/bin/gatdaem1d.dll")
		message(STATUS "Copying shared library gatdaem1d.dll to python directory")
		configure_file("${CMAKE_INSTALL_PREFIX}/bin/gatdaem1d.dll" "${CMAKE_INSTALL_PREFIX}/python/gatdaem1d/gatdaem1d.dll" COPYONLY)
		message(STATUS "Copying shared library gatdaem1d.dll to matlab directory as gatdaem1d.mexw64")
		configure_file("${CMAKE_INSTALL_PREFIX}/bin/gatdaem1d.dll" "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.mexw64" COPYONLY)
	endif()
endif()




