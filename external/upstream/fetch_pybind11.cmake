#set(PYBIND11_PYTHON_VERSION 3.6)
#set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")

find_package(pybind11 2.6 CONFIG QUIET)
if(pybind11_FOUND)
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR} (found version ${pybind11_VERSION})")
else()
  message(STATUS "Suitable pybind11 could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(pybind11
    QUIET
    URL
      https://github.com/pybind/pybind11/archive/v2.6.0.tar.gz
    )
  #set(PYBIND11_PYTHON_VERSION ${PYBIND11_PYTHON_VERSION})
  set(PYBIND11_TEST OFF CACHE BOOL "")
  FetchContent_MakeAvailable(pybind11)
endif()
