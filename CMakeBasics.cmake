find_package(PythonInterp)
function(instantiate _inst_var _source_var)
  set(_dims 1)
  if(ENABLE_2D)
    list(APPEND _dims 2)
  endif()
   if(ENABLE_3D)
    list(APPEND _dims 3)
  endif()

  set(_tmp_sources )
  foreach(_inst ${${_inst_var}})
    if(ENABLE_INSTANTIATIONS)
      foreach(_dim ${_dims})
        set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
        set(_inst_file ${_inst}${_dim}d.cc)
        #Run python script to generate instantiation file
        execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/helpers/InstantiationGenerator.py 
                                                                          ${_inst_py} ${_inst_file} ${_dim}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                        )
        list(APPEND _tmp_sources ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
      endforeach()
    else()
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_inst}.cc")
        list(APPEND _tmp_sources ${_inst}.cc)
      else()
        message("WARNING : Could not find ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}.cc")
      endif()      
    endif()      
  endforeach()
  set(${_source_var} ${${_source_var}} ${_tmp_sources} PARENT_SCOPE)
endfunction()