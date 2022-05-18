# Export compile commands for each file to JSON
# This is useful for static analysis tools and linters
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED COMPONENTS CXX)

find_package(OpenMP 4.5 REQUIRED COMPONENTS CXX)

message(STATUS "Configuring a ${CMAKE_BUILD_TYPE} build")
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "C++ compiler flags")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
message(STATUS "   Build-type-specific    : ${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}")
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(STATUS "   Project defaults       : ${VLX_CXX_FLAGS}")
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")
message(STATUS "   OpenMP parallelization : ${OpenMP_CXX_FLAGS}")

# transform VLX_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
string(REPLACE " " ";" VLX_CXX_FLAGS ${VLX_CXX_FLAGS})
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS ${EXTRA_CXXFLAGS})
endif()

# figure out where to put the Python module
if(NOT DEFINED PYMOD_INSTALL_FULLDIR)
  if(NOT WIN32)
    # install Python module under CMAKE_INSTALL_LIBDIR
    # if that is "lib64", the use just "lib"
    set(_lib "${CMAKE_INSTALL_LIBDIR}")
    if(CMAKE_INSTALL_LIBDIR STREQUAL "lib64")
      set(_lib "lib")
    endif()
    set(PYMOD_INSTALL_FULLDIR
          "${_lib}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/veloxchem"
        CACHE STRING "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed" FORCE)
  else()
    set(PYMOD_INSTALL_FULLDIR
          "Lib/site-packages/veloxchem"
        CACHE STRING "Location under CMAKE_INSTALL_PREFIX to which Python modules are installed" FORCE)
  endif()
endif()
message(STATUS "Setting PYMOD_INSTALL_FULLDIR: ${PYMOD_INSTALL_FULLDIR}")
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

add_subdirectory(src)

# handle folder with basis sets
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/basis)
# we glob the basis set files in basis and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
  GLOB
    _vlx_basis_sets
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/basis/*
  )
# 1. symlink under the build tree
foreach(_basis IN LISTS _vlx_basis_sets)
 get_filename_component(__basis ${_basis} NAME)
 file(
   CREATE_LINK
     ${_basis}
     ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/basis/${__basis}
   COPY_ON_ERROR
   SYMBOLIC
   )
endforeach()
# 2. install rules for basis sets folder
install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )

enable_testing()
include(CTest)
# This must come last!!
add_subdirectory(unit_tests)
add_subdirectory(python_tests)
