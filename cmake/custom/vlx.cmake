# Export compile commands for each file to JSON
# This is useful for static analysis tools and linters
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# OpenMP is found in the top-level CMakeLists.txt, before the compiler flags
# module, because the vector math library probe links against
# OpenMP::OpenMP_CXX.

# The non-SIMD three-center electron repulsion drivers (and the momentum
# transfer relations they are built on) are superseded by the SIMD
# implementation and are not used anywhere outside the Python bindings.
# They are therefore not compiled unless explicitly requested.
option_with_print(VLX_USE_OLD_RI "Build the deprecated non-SIMD three-center RI integral drivers" OFF)

message(STATUS "Configuring a ${CMAKE_BUILD_TYPE} build")
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "C++ compiler flags")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
message(STATUS "   Build-type-specific    : ${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}")
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(STATUS "   Project defaults       : ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${VLX_CXX_FLAGS}")
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")
message(STATUS "   OpenMP parallelization : ${OpenMP_CXX_FLAGS}")

# transform VLX_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
string(REPLACE " " ";" VLX_CXX_FLAGS "${VLX_CXX_FLAGS}")
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS "${EXTRA_CXXFLAGS}")
endif()

# figure out where to put the Python module
if(NOT DEFINED PYMOD_INSTALL_FULLDIR)
  if(NOT WIN32)
    set(PYMOD_INSTALL_FULLDIR
          "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/veloxchem"
        CACHE STRING
          "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed"
        FORCE
      )
  else()
    set(PYMOD_INSTALL_FULLDIR
          "Lib/site-packages/veloxchem"
        CACHE STRING
          "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed"
        FORCE
      )
  endif()
endif()
message(STATUS "Setting PYMOD_INSTALL_FULLDIR: ${PYMOD_INSTALL_FULLDIR}")
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

add_subdirectory(src)

# handle folder with data files
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/database
        DESTINATION ${PYMOD_INSTALL_FULLDIR}
        FILES_MATCHING PATTERN "*")

# handle folder with basis files
install(DIRECTORY   ${PROJECT_SOURCE_DIR}/basis
        DESTINATION ${PYMOD_INSTALL_FULLDIR}
        FILES_MATCHING PATTERN "*")
