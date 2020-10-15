# find mpi4py
include(${PROJECT_SOURCE_DIR}/cmake/custom/FindPythonModule.cmake)

find_python_module(mpi4py REQUIRED)

# we also need the include directories for mpi4py
if(mpi4py_FOUND)
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" "-c"
      "import mpi4py as m; print(m.__version__); print(m.get_include());"
    RESULT_VARIABLE
      _mpi4py_SEARCH_SUCCESS
    OUTPUT_VARIABLE
      _mpi4py_VALUES
    ERROR_VARIABLE
      _mpi4py_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  # Convert the process output into a list
  string(REGEX REPLACE ";" "\\\\;" _mpi4py_VALUES ${_mpi4py_VALUES})
  string(REGEX REPLACE "\n" ";" _mpi4py_VALUES ${_mpi4py_VALUES})
  list(GET _mpi4py_VALUES 0 mpi4py_VERSION)
  list(GET _mpi4py_VALUES 1 mpi4py_INCLUDE_DIRS)

  # Make sure all directory separators are '/'
  string(REGEX REPLACE "\\\\" "/" mpi4py_INCLUDE_DIRS ${mpi4py_INCLUDE_DIRS})

  # Get the major and minor version numbers
  string(REGEX REPLACE "\\." ";" _mpi4py_VERSION_LIST ${mpi4py_VERSION})
  list(GET _mpi4py_VERSION_LIST 0 mpi4py_VERSION_MAJOR)
  list(GET _mpi4py_VERSION_LIST 1 mpi4py_VERSION_MINOR)
  list(GET _mpi4py_VERSION_LIST 2 mpi4py_VERSION_PATCH)
  string(REGEX MATCH "[0-9]*" mpi4py_VERSION_PATCH ${mpi4py_VERSION_PATCH})
  math(EXPR mpi4py_VERSION_DECIMAL
      "(${mpi4py_VERSION_MAJOR} * 10000) + (${mpi4py_VERSION_MINOR} * 100) + ${mpi4py_VERSION_PATCH}")
endif()

# create Python module with pybind11
# * use pybind11 header files with -isystem
# * explicitly disable LTO, even for release builds
pybind11_add_module(veloxchemlib
  MODULE SYSTEM NO_EXTRAS
    ${CMAKE_CURRENT_LIST_DIR}/ExportDFT.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportGeneral.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportMath.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportMolData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportOneInts.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportOrbData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportResponse.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportTwoInts.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportVeloxChemLib.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportVisualization.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExportXTB.cpp
  )

target_precompile_headers(veloxchemlib
  PRIVATE
    <mpi4py/mpi4py.h>
    <pybind11/numpy.h>
    <pybind11/operators.h>
    <pybind11/pybind11.h>
    <pybind11/stl.h>
  )

target_include_directories(veloxchemlib
  SYSTEM PRIVATE
    ${mpi4py_INCLUDE_DIRS}
  )

target_link_libraries(veloxchemlib
  PRIVATE
    vlxobjs
  )

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(VLX_RPATH "@loader_path/${_rel}")
else()
  set(VLX_RPATH "\$ORIGIN/${_rel}")
endif()

set_target_properties(veloxchemlib
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${VLX_RPATH}${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
    RESOURCE "${_veloxchemlib_pys}"
  )

install(
  TARGETS
    veloxchemlib
  ARCHIVE
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  LIBRARY
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  RUNTIME
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  RESOURCE
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    COMPONENT lib
  )
