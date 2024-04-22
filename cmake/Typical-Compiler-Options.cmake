# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Set the compile flags
if(MSVC)
	add_compile_options(/O2)
	add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
	add_compile_definitions(_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
else()
	add_compile_options(-O3)
	add_compile_options(-Wno-unused-but-set-variable)
	add_compile_options(-Wno-sign-compare)
	add_compile_options(-Wno-format-security)
	add_compile_options(-Wno-tautological-constant-compare)
	link_libraries(-lstdc++fs)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
	add_compile_options(-Wno-unused-result)
	add_compile_options(-Wno-date-time)
	add_compile_options(-Wno-error=date-time)
	#On GCC, even with -Wno-date-time, still get warings of the form: warning: macro "__DATE__" might prevent reproducible builds [-Wdate-time]
endif()
