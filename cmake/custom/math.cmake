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
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ``LA_FOUND``
#   library implementing the BLAS and LAPACK interfaces is found
# ``LA_LINKER_FLAGS``
#   uncached list of required linker flags (excluding ``-l`` and ``-L``).
# ``LA_LIBRARIES``
#   uncached list of libraries (using full path name) to link against
#   to use BLAS and LAPACK


option_with_default(VLX_LA_VENDOR "Linear algebra library vendor" "MKL")
list(APPEND _valid_vendors "MKL" "OpenBLAS")
if(DEFINED VLX_LA_VENDOR AND NOT VLX_LA_VENDOR IN_LIST _valid_vendors)
  message(STATUS "${VLX_LA_VENDOR} not a valid alignment, resetting to default value MKL")
  set(VLX_LA_VENDOR "MKL" CACHE STRING "Linear algebra vendor" FORCE)
endif()

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
      "${_mkl_path_suffixes}"
    )
  if(CXX_COMPILER_ID MATCHES GNU OR CXX_COMPILER_ID MATCHES AppleClang)
    set(_la_compiler_flags "-m64")
  endif()
else()
  # locate cblas.h and lapacke.h
  find_path(_la_include_dirs
    NAMES
      cblas.h
      lapacke.h
    HINTS
      ${CMAKE_INSTALL_PREFIX}/include
    )
endif()

if(NOT _la_include_dirs)
  message(FATAL_ERROR "Could not find header files for ${VLX_LA_VENDOR} linear algebra backend")
endif()

find_package(LAPACK REQUIRED)

# create target
# adapted from: https://github.com/Kitware/CMake/blob/master/Modules/FindLAPACK.cmake#L588-L603
if(BLAS_FOUND AND LAPACK_FOUND AND NOT TARGET Math::LA)
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
  endif()
  unset(_la_compile_definitions)
  unset(_la_compiler_flags)
  unset(_la_include_dirs)
  unset(_la_libs)
endif()
