find_package(PythonInterp)

function(instantiate _inst_var _source_var)
  set(_dims)
  if(ENABLE_2D)
    list(APPEND _dims 2)
  endif()
   if(ENABLE_3D)
    list(APPEND _dims 3)
  endif()

  set(_tmp_sources)
  foreach(_inst ${${_inst_var}})
    if(ENABLE_INSTANTIATIONS)
      foreach(_dim ${_dims})
	set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
	set(_inst_file ${_inst}Inst.cc.${_dim}.cc)
	#Run python script to generate instantiation file
	execute_process(COMMAND ${PYTHON_EXECURTABLE} ${PROJECT_NAME}/src/helpers/InstantiationGenerator.py 
									  ${_inst_py} ${_inst_file} ${_dim}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			)
	list(APPEND _tmp_sources ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
      endforeach()
    else()
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_file}.cc)
	list(APPEND _tmp_source ${_inst}.cc)
      endif()
    endif()      
  endforeach()
  set(${_source_var} ${${_source_var}} ${_tmp_source} PARENT_SCOPE)
endfunction()
