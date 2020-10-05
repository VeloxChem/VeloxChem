set(PYBIND11_PYTHON_VERSION 3.6)
set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")

find_package(pybind11 2.5 CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(pybind11_sources
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.5.0.tar.gz
    )

  FetchContent_GetProperties(pybind11_sources)

  set(PYBIND11_PYTHON_VERSION ${PYBIND11_PYTHON_VERSION})
  set(PYBIND11_TEST OFF CACHE BOOL "")
  # FIXME Needed?
  #set(PYMOD_INSTALL_FULLDIR ${PYMOD_INSTALL_FULLDIR})

  if(NOT pybind11_sources_POPULATED)
    FetchContent_Populate(pybind11_sources)

    add_subdirectory(
      ${pybind11_sources_SOURCE_DIR}
      ${pybind11_sources_BINARY_DIR}
      )
  endif()
endif()
