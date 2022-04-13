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
#   The list of valid vendors in this module is the same as for the
#   official ``FindBLAS.cmake`` module, plus ``Cray``.
#   Note that this module reimplements the detection of ``FLAME`` to be able to
#   discover its multithreaded version.
#   Defaults to ``Generic``.
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``Math::LA``
#   The libraries to use as linear algebra backend (BLAS and LAPACK), if found.
#

option_with_default(VLX_LA_VENDOR "Linear algebra library vendor" "Generic")

# we use the standard CMake FindLapack module
# it will search for BLAS too
# we set up the call by defining the BLA_VENDOR variable
set(BLA_VENDOR ${VLX_LA_VENDOR})

# compile definitions for the Math::LA target
list(APPEND _la_compile_definitions
  "HAS_${VLX_LA_VENDOR}"
  "ENABLE_${VLX_LA_VENDOR}"
  "VLX_HAS_${VLX_LA_VENDOR}"
  )

# compiler flags for the Math::LA target
set(_la_compiler_flags)

# include directories for the Math::LA target
set(_la_include_dirs)

# libraries for the Math::LA target
set(_la_libraries)


if(VLX_LA_VENDOR MATCHES "MKL")
  # MKL is known under a different name to CMake
  set(BLA_VENDOR "Intel10_64_dyn")
  find_package(LAPACK REQUIRED)
  list(APPEND _la_libraries ${LAPACK_LIBRARIES})

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
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${_mkl_root}
    PATH_SUFFIXES
      ${_mkl_path_suffixes}
    )

  if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(_la_compiler_flags "-m64")
  endif()
elseif(VLX_LA_VENDOR MATCHES "FLAME")
  # detection of BLIS/FLAME in the standard CMake module cannot find the
  # multithreaded version, hence the code here.
  include(CheckFunctionExists)
  include(FindPackageHandleStandardArgs)

  list(APPEND _extaddlibdir ENV LD_LIBRARY_PATH "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  # we would prefer the multithreaded version, if available
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

  find_package_handle_standard_args(BLAS
    REQUIRED_VARS
      _la_blis_library
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

  find_package_handle_standard_args(LAPACK
    REQUIRED_VARS
      _la_flame_library
      _la_flame_library_WORKS
    )

  list(APPEND _la_libraries ${_la_flame_library} ${_la_blis_library})
  list(APPEND _la_include_dirs
    ${PROJECT_SOURCE_DIR}/external/upstream/cblas
    ${PROJECT_SOURCE_DIR}/external/upstream/lapacke
    )

  unset(_la_blis_library)
  unset(_la_flame_library)
elseif(VLX_LA_VENDOR MATCHES "Cray")
  # compiling and linking with Cray libsci is taken care of automatically by the
  # Cray compiler wrappers. Do nothing!
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
else()
  find_package(LAPACK REQUIRED)
  list(APPEND _la_libraries ${LAPACK_LIBRARIES})
  list(APPEND _la_include_dirs
    ${PROJECT_SOURCE_DIR}/external/upstream/cblas
    ${PROJECT_SOURCE_DIR}/external/upstream/lapacke
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

  target_compile_definitions(Math::LA
    INTERFACE
      "${_la_compile_definitions}"
    )

  target_compile_options(Math::LA
    INTERFACE
      "${_la_compiler_flags}"
    )

  target_include_directories(Math::LA
    SYSTEM
    INTERFACE
      "${_la_include_dirs}"
    )

  target_link_libraries(Math::LA
    INTERFACE
      "${_la_libraries}"
  )

  unset(_la_compile_definitions)
  unset(_la_compiler_flags)
  unset(_la_include_dirs)
  unset(_la_libraries)
endif()

# print out some info on the Math::LA target
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
