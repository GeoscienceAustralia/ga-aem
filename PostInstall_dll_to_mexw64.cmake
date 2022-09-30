message(STATUS "Running cmake post install script")
message(STATUS "Renaming matlab .dll to .mexw64")
file(RENAME "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.dll" "${CMAKE_INSTALL_PREFIX}/matlab/bin/gatdaem1d.mexw64")

