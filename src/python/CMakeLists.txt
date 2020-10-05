target_sources(veloxchemlib
  PRIVATE
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
  )

# find mpi4py, this is always needed, even when building without MPI
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

target_include_directories(veloxchemlib
  SYSTEM PRIVATE
    ${mpi4py_INCLUDE_DIRS}
  )
