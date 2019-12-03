#-----------------------------------------------------------------------------------
# - A Python-wrapping module for CMake using PYB11Generator/pybind11
#
# Modified from the UseSWIG module from the cmake2.8 release.
#
# Defines the following macros:
#   PYB11_GENERATE_BINDINGS(module_name module_list)
#     - Generates the Python bindings for each module in the list
#
# Internal macros to do the work:
#   PYB11_REGISTER_MODULE(name)
#     - Builds list of Python source code and generated C++ code
# 
#
# Variables that must be set before calling PYB11_GENERATE_BINDINGS:
#   PYB11_DIR
#     - Location of the python source code to build the bindings
#   PYTHON_EXECUTABLE
#     - Python executable
#   PYTHON_LIB_DIR
#     - Python lib (typically obtained from the executable root)
#
# Optional variables
#   PYB11_ADDITIONAL_ARGS
#     - Additional arguments to follow the command in PYB11_GENERATE_BINDINGS
#
# To get the names of the generated source
# use: ${PYB11_GENERATED_SOURCE}
#-----------------------------------------------------------------------------------

# Initialize the list of Python source code
set(PYB11_SOURCE)

# Initialize the list of C++ code generated by PYB11Generator
set(PYB11_GENERATED_SOURCE)

#
# INTERNAL: Build up the lists of source files and generated files to track
# outputs and dependencies for the build
#
macro(PYB11_REGISTER_MODULE _module)
  if (EXISTS "${PYB11_DIR}/${_module}MOD.py")
    set(_source "${PYB11_DIR}/${_module}MOD.py")
  else()
    if (EXISTS "${PROJECT_BINARY_DIR}/src/PYB11/${_module}MOD.py")
      set(_source "${PROJECT_BINARY_DIR}/src/PYB11/${_module}MOD.py")
    else()
      message(FATAL_ERROR "Could not find source file ${_module}MOD.py")
    endif()
  endif()
  list(APPEND PYB11_SOURCE ${_source})
  list(APPEND PYB11_GENERATED_SOURCE "${PYB11_MODULE_NAME}_${_module}.cc")
endmacro()

#
# Runs "python polytopeMOD.py" to generate the bindings.
# This is where the PYB11Generator module is actually used. Calls the previous
# internal macros to generate the source file lists
#
macro(PYB11_GENERATE_BINDINGS)
  set(PYB11_GENERATED_SOURCE)
  list(APPEND PYB11_SOURCE "${PYB11_MODULE_NAME}MOD.py")
  list(APPEND PYB11_GENERATED_SOURCE "${PYB11_MODULE_NAME}.cc")

  set(PYTHON_ENV 
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/polytope:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Distributed:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/OpenMP:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/CXXTypes:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Geometry:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/PolyClipper:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Silo:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/DataOutput:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/NodeList:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Field:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Kernel:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Neighbor:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Material:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/FileIO:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/DataBase:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Boundary:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Physics:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Hydro:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ExternalForce:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Gravity:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Integrator:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Utilities:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/NodeGenerators:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/FieldOperations:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SPH:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/CRKSPH:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ArtificialViscosity:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SVPH:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Mesh:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Damage:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SolidMaterial:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Strength:"
    "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ArtificialConduction:"
    "${PROJECT_SOURCE_DIR}/src/SimulationControl:"
    )
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})

  add_custom_command(
    OUTPUT Spheral${PYB11_GENERATED_SOURCE} ${PYB11_GENERATED_HEADER}
    COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
    ${PYTHON_EXECUTABLE} -c
    'from PYB11Generator import * \; 
    import ${PYB11_MODULE_NAME}MOD \;
    PYB11generateModule(${PYB11_MODULE_NAME}MOD, \"Spheral${PYB11_MODULE_NAME}\") '
    DEPENDS ${PYB11_SOURCE}
    )
endmacro()