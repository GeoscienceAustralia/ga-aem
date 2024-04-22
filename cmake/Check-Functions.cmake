function(check_existance invar rstatus)
	set(${rstatus} TRUE PARENT_SCOPE)
	foreach(item ${${invar}})
		cmake_path(SET item NORMALIZE "${item}")
		if (EXISTS "${item}")
			#message(STATUS "${dir} exists")
		else()
			#message(STATUS "${dir} does not exist")
			set(${rstatus} FALSE PARENT_SCOPE)
		endif()
	endforeach()
endfunction()


