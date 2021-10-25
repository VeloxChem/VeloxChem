# make sure site-packages is on path
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE
    _tmp_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_QUIET
  ERROR_QUIET
  )

find_package(pybind11 2.6 CONFIG QUIET
  HINTS
    ${_tmp_dir}
  )

if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(pybind11
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.8.0.tar.gz
    )
  #set(PYBIND11_PYTHON_VERSION ${PYBIND11_PYTHON_VERSION})
  set(PYBIND11_TEST OFF CACHE BOOL "")
  FetchContent_MakeAvailable(pybind11)
endif()

unset(_tmp_dir)
