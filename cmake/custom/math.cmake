#.rst:
#
# Find linear algebra backend libraries (BLAS and LAPACK)
#
# Input Variables
# ^^^^^^^^^^^^^^^
#
# The following variables may be set to influence this module's behavior:
#
# ``VLX_LA_VENDOR``
#   Check if libraries are available for the given vendor.
#   List of vendors valid in this module:
#
#   * ``OpenBLAS``
#   * ``MKL``
#   * ``Cray``
#   * ``FLAME``
#
#   Defaults to ``MKL``.
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``Math::LA``
#   The libraries to use as linear algebra backend (BLAS and LAPACK), if found.
#

option_with_default(VLX_LA_VENDOR "Linear algebra library vendor" "MKL")
list(APPEND _valid_vendors "MKL" "OpenBLAS" "Cray" "FLAME")

if(DEFINED VLX_LA_VENDOR AND NOT VLX_LA_VENDOR IN_LIST _valid_vendors)
  message(STATUS "${VLX_LA_VENDOR} not a valid vendor, resetting to default value MKL")
  set(VLX_LA_VENDOR "MKL" CACHE STRING "Linear algebra vendor" FORCE)
endif()

include(CheckFunctionExists)
include(FindPackageHandleStandardArgs)

# we use the standard CMake FindLapack module
# it will search for BLAS too
# we set up the call by defining the BLA_VENDOR variable
set(BLA_VENDOR ${VLX_LA_VENDOR})
list(APPEND _la_compile_definitions "HAS_${VLX_LA_VENDOR}" "ENABLE_${VLX_LA_VENDOR}" "VLX_HAS_${VLX_LA_VENDOR}")
set(_la_compiler_flags "")

if(VLX_LA_VENDOR MATCHES "MKL")
  set(BLA_VENDOR "Intel10_64_dyn")
  # locate mkl.h
  # MKL uses a multitude of partially platform-specific subdirectories:
  if(WIN32)
    set(LAPACK_mkl_OS_NAME "win")
  elseif(APPLE)
    set(LAPACK_mkl_OS_NAME "mac")
  else()
    set(LAPACK_mkl_OS_NAME "lin")
  endif()
  if(DEFINED ENV{MKLROOT})
    file(TO_CMAKE_PATH "$ENV{MKLROOT}" _mkl_root)
    # If MKLROOT points to the subdirectory 'mkl', use the parent directory instead
    # so we can better detect other relevant libraries in 'compiler' or 'tbb':
    get_filename_component(_mkl_root_last_dir "${_mkl_root}" NAME)
    if(_mkl_root_last_dir STREQUAL "mkl")
      get_filename_component(_mkl_root "${_mkl_root}" DIRECTORY)
    endif()
  endif()
  list(APPEND _mkl_path_suffixes
    "include"
    "compiler/include"
    "compiler/include/intel64_${LAPACK_mkl_OS_NAME}"
    "compiler/include/intel64"
    "mkl/include" "mkl/include/intel64_${LAPACK_mkl_OS_NAME}"
    "mkl/include/intel64"
    "include/intel64_${LAPACK_mkl_OS_NAME}"
    )

  find_path(_la_include_dirs
    NAMES
      mkl.h
    HINTS
      ${CMAKE_INSTALL_PREFIX}/include
    PATHS
      ${_mkl_root}
    PATH_SUFFIXES
      ${_mkl_path_suffixes}
    )
  if(CXX_COMPILER_ID MATCHES GNU OR CXX_COMPILER_ID MATCHES AppleClang)
    set(_la_compiler_flags "-m64")
  endif()

  find_package(LAPACK REQUIRED)
elseif(VLX_LA_VENDOR MATCHES "FLAME")
  # CMake detection of BLIS and FLAME, since it can't handle their multithreaded version
  # for BLIS, we would prefere the multithreaded version, if available
  list(APPEND _extaddlibdir ENV LD_LIBRARY_PATH "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  find_library(_la_blis_library
    NAMES
      "blis-mt"
      "blis"
    NAMES_PER_DIR
    PATHS
      ${_extaddlibdir}
    )
  set(CMAKE_REQUIRED_LIBRARIES ${_la_blis_library})
  check_function_exists("sgemm_" _la_blis_library_WORKS)
  set(CMAKE_REQUIRED_LIBRARIES)

  # locate cblas.h
  # get the directory containing the library
  get_filename_component(_hint ${_la_blis_library} DIRECTORY)
  # replace "lib" (or "lib64") by "include"
  string(REGEX REPLACE
    "(lib|lib64)"
    "include/blis"
    _hint
    ${_hint}
    )

  find_path(_la_blis_include_dirs
    NAMES
      cblas.h
    HINTS
      ${_hint}
      ${CMAKE_INSTALL_PREFIX}/include
    )

  find_package_handle_standard_args(BLAS
    REQUIRED_VARS
      _la_blis_library
      _la_blis_include_dirs
      _la_blis_library_WORKS
    )

  find_library(_la_flame_library
    NAMES
      "flame"
    NAMES_PER_DIR
    PATHS
      ${_extaddlibdir}
    )
  set(CMAKE_REQUIRED_LIBRARIES ${_la_flame_library})
  check_function_exists("cheev_" _la_flame_library_WORKS)
  set(CMAKE_REQUIRED_LIBRARIES)

  # locate lapacke.h
  # get the directory containing the library
  get_filename_component(_hint ${_la_flame_library} DIRECTORY)
  # replace "lib" (or "lib64") by "include"
  string(REGEX REPLACE
    "(lib|lib64)$"
    "include"
    _hint
    ${_hint}
    )

  find_path(_la_flame_include_dirs
    NAMES
      lapacke.h
    HINTS
      ${_hint}
      ${CMAKE_INSTALL_PREFIX}/include
    )

  find_package_handle_standard_args(LAPACK
    REQUIRED_VARS
      _la_flame_library
      _la_flame_include_dirs
      _la_flame_library_WORKS
    )

  list(APPEND LAPACK_LIBRARIES ${_la_flame_library} ${_la_blis_library})
  list(APPEND _la_include_dirs ${_la_flame_include_dirs} ${_la_blis_include_dirs})

  unset(_la_blis_library)
  unset(_la_flame_library)
  unset(_la_blis_include_dirs)
  unset(_la_flame_include_dirs)
elseif(VLX_LA_VENDOR MATCHES "Cray")
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
else()
  find_package(LAPACK REQUIRED)

  # locate cblas.h and lapacke.h
  # to create a hint for locating the headers, we extract the root of the
  # installation from the last element in LAPACK_LIBRARIES
  list(POP_BACK LAPACK_LIBRARIES _hint)
  # get the directory containing the library
  get_filename_component(_hint ${_hint} DIRECTORY)
  # replace "lib" (or "lib64") by "include"
  string(REGEX REPLACE 
    "(lib|lib64)"
    "include"
    _hint 
    ${_hint}
    )

  find_path(_la_include_dirs
    NAMES
      cblas.h
      lapacke.h
    HINTS
      ${_hint}
      ${CMAKE_INSTALL_PREFIX}/include
    )
endif()

if(NOT _la_include_dirs AND NOT VLX_LA_VENDOR MATCHES "Cray")
  message(FATAL_ERROR "Could not find header files for ${VLX_LA_VENDOR} linear algebra backend")
endif()

# create target
# adapted from: https://github.com/Kitware/CMake/blob/master/Modules/FindLAPACK.cmake#L588-L603
if(BLAS_FOUND AND LAPACK_FOUND AND NOT TARGET Math::LA)
  message(STATUS "Using ${VLX_LA_VENDOR} as linear algebra backend")
  add_library(Math::LA INTERFACE IMPORTED)
  set(_la_libs "${LAPACK_LIBRARIES}")

  if(_la_libs)
    set_target_properties(Math::LA
      PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "${_la_compile_definitions}"
        INTERFACE_COMPILE_OPTIONS "${_la_compiler_flags}"
        INTERFACE_INCLUDE_DIRECTORIES "${_la_include_dirs}"
        INTERFACE_LINK_LIBRARIES "${_la_libs}"
    )
    include(CMakePrintHelpers)
    cmake_print_properties(
      TARGETS 
        Math::LA
      PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS
        INTERFACE_COMPILE_OPTIONS
        INTERFACE_INCLUDE_DIRECTORIES
        INTERFACE_LINK_LIBRARIES
      )
  endif()
  unset(_la_compile_definitions)
  unset(_la_compiler_flags)
  unset(_la_include_dirs)
  unset(_la_libs)
endif()
