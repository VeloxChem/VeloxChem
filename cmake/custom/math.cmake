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
#
# Notes
# ^^^^^
#
# We offer configuration shortcuts for MKL and OpenBLAS. Setting the ``MKLROOT``
# and ``OPENBLASROOT`` environment variables will be enough to configure the
# code to link against and use MKL and OpenBLAS, respectively.
#
# This module reimplements the detection of ``FLAME`` to be able to also
# discover the multithreaded versions of both the BLIS and FLAME libraries.
#
# The ``Generic`` option is available to link against the reference Netlib
# implementation of BLAS and LAPACK.
# Essentially, it is meant to be used for building the Conda package and switch
# the BLAS/LAPACK implementation to use at **installation time**.  It should
# **never** be used in production!
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``Math::LA``
#   The libraries to use as linear algebra backend (BLAS and LAPACK), if found.

if(NOT DEFINED VLX_LA_VENDOR OR "${VLX_LA_VENDOR}" STREQUAL "")
  if(DEFINED ENV{MKLROOT})
    set(VLX_LA_VENDOR "MKL")
    file(TO_CMAKE_PATH "$ENV{MKLROOT}" _root)
    set(LAPACK_ROOT ${_root})
  elseif(DEFINED ENV{OPENBLASROOT})
    set(VLX_LA_VENDOR "OpenBLAS")
    file(TO_CMAKE_PATH "$ENV{OPENBLASROOT}" _root)
    set(LAPACK_ROOT ${_root})
  else()
    message(FATAL_ERROR "VLX_LA_VENDOR undefined! Please, let me know with linear algebra backend you'd like to use.")
  endif()
endif()

print_option(VLX_LA_VENDOR "${VLX_LA_VENDOR}")
set(VLX_LA_VENDOR "${VLX_LA_VENDOR}" CACHE STRING "Linear algebra library vendor" FORCE)

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
set(_la_compiler_flags "")

# include directories for the Math::LA target
set(_la_include_dirs "")

# whether to check that CBLAS and LAPACKE work
set(_check_cblas_and_lapacke TRUE)

if(VLX_LA_VENDOR STREQUAL "MKL")
  # MKL is known under a different name to CMake
  set(BLA_VENDOR "Intel10_64_dyn")
  find_package(LAPACK REQUIRED QUIET)

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

  find_path(_mkl_include_dirs
    NAMES
      mkl.h
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${_mkl_root}
    PATH_SUFFIXES
      ${_mkl_path_suffixes}
    )
  # clear generic _la_include_dirs and use MKL own
  set(_la_include_dirs "${_mkl_include_dirs}")

  if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(_la_compiler_flags "-m64")
  endif()

  set(_check_cblas_and_lapacke FALSE)
elseif(VLX_LA_VENDOR STREQUAL "Apple" OR VLX_LA_VENDOR STREQUAL "NAS")
  set(CMAKE_THREAD_LIBS_INIT "-lpthread")
  set(CMAKE_HAVE_THREADS_LIBRARY 1)
  set(CMAKE_USE_WIN32_THREADS_INIT 0)
  set(CMAKE_USE_PTHREADS_INIT 1)
  set(THREADS_PREFER_PTHREAD_FLAG ON)

  find_package(LAPACK REQUIRED QUIET)

  set(_check_cblas_and_lapacke FALSE)
elseif(VLX_LA_VENDOR STREQUAL "Cray")
  # compiling and linking with Cray libsci is taken care of automatically by the
  # Cray compiler wrappers. Do nothing!
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
  # clear generic _la_include_dirs and let the Cray compiler figure it out
  set(_la_include_dirs "")
  set(LAPACK_LIBRARIES "")

  set(_check_cblas_and_lapacke FALSE)
elseif(VLX_LA_VENDOR STREQUAL "FLAME")
  # detection of BLIS/FLAME in the standard CMake module cannot find the
  # multithreaded version, hence the code here.
  include(FindPackageHandleStandardArgs)
  include(CheckCXXSymbolExists)

  # check_* invocations will not print any output
  set(CMAKE_REQUIRED_QUIET ON)

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
  check_cxx_symbol_exists("cblas_sgemm" "cblas.h" _la_blis_library_works)
  set(CMAKE_REQUIRED_LIBRARIES)

  find_package_handle_standard_args(BLAS
    REQUIRED_VARS
      _la_blis_library
      _la_blis_library_works
    )

  find_library(_la_flame_library
    NAMES
      "flame"
    NAMES_PER_DIR
    PATHS
      ${_extaddlibdir}
    )
  set(CMAKE_REQUIRED_LIBRARIES ${_la_flame_library})
  check_cxx_symbol_exists("LAPACKE_cheev" "lapacke.h" _la_flame_library_works)
  set(CMAKE_REQUIRED_LIBRARIES)

  find_package_handle_standard_args(LAPACK
    REQUIRED_VARS
      _la_flame_library
      _la_flame_library_works
    )

  list(APPEND LAPACK_LIBRARIES ${_la_flame_library} ${_la_blis_library})

  unset(_la_blis_library)
  unset(_la_blis_library_works)
  unset(_la_flame_library)
  unset(_la_flame_library_works)

  set(CMAKE_REQUIRED_QUIET OFF)

  set(_check_cblas_and_lapacke FALSE)
elseif(VLX_LA_VENDOR STREQUAL "Generic")
  find_package(LAPACK REQUIRED QUIET)

  set(_check_cblas_and_lapacke TRUE)
else()
  find_package(LAPACK REQUIRED QUIET)

  set(_check_cblas_and_lapacke TRUE)
endif()

# check_* invocations will not print any output
set(CMAKE_REQUIRED_QUIET ON)

# check that CBLAS works
if(_check_cblas_and_lapacke)
  include(${CMAKE_CURRENT_LIST_DIR}/cblas-lapacke.cmake)

  # get directory containing cblas.h
  # we are assuming this is on an already-known system path
  find_path(_include_dirs
    NAMES
      cblas.h
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${BLAS_ROOT}/include
      ${LAPACK_ROOT}/include
    )
  list(APPEND _la_include_dirs ${_include_dirs})

  # is CBLAS bundled with the BLAS_LIBRARIES we found?
  set(_blas_libs)
  foreach(_lib IN LISTS BLAS_LIBRARIES)
    get_filename_component(__lib ${_lib} NAME_WE)
    list(APPEND _blas_libs ${__lib})
  endforeach()
  list(TRANSFORM _blas_libs REPLACE "^lib" "")

  check_cblas_lapacke(_la_cblas_bundled
    CBLAS
    cblas_sgemm
    ${_include_dirs}
    "cblas.h"
    ""
    "${_blas_libs}"
    ""
    "${BLAS_ROOT};${LAPACK_ROOT}"
    "lib"
    )

  if(NOT _la_cblas_bundled)
    # well, it's not bundled, let's give it another shot by also linking the
    # generic libcblas.so
    check_cblas_lapacke(_la_cblas_library
      CBLAS
      cblas_sgemm
      ${_include_dirs}
      "cblas.h"
      ""
      "cblas;${_blas_libs}"
      ""
      "${BLAS_ROOT};${LAPACK_ROOT}"
      "lib"
      )

    if(_la_cblas_library)
      list(PREPEND BLAS_LIBRARIES ${_la_cblas_library})
      list(INSERT LAPACK_LIBRARIES 1 ${_la_cblas_library})
      unset(_la_cblas_library)
      unset(_include_dirs)
    else()
      message(FATAL_ERROR "Cannot find a working CBLAS library to link to!")
    endif()
  endif()

  # get directory containing lapacke.h
  # we are assuming this is on an already-known system path
  find_path(_include_dirs
    NAMES
      lapacke.h
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${BLAS_ROOT}/include
      ${LAPACK_ROOT}/include
    )
  list(APPEND _la_include_dirs ${_include_dirs})

  # is LAPACKE bundled with the LAPACK_LIBRARIES we found?
  set(_lapack_libs)
  foreach(_lib IN LISTS LAPACK_LIBRARIES)
    get_filename_component(__lib ${_lib} NAME_WE)
    list(APPEND _lapack_libs ${__lib})
  endforeach()
  list(TRANSFORM _lapack_libs REPLACE "^lib" "")

  check_cblas_lapacke(_la_lapacke_bundled
    LAPACKE
    LAPACKE_cheev
    ${_include_dirs}
    "lapacke.h"
    ""
    "${_lapack_libs}"
    ""
    "${BLAS_ROOT};${LAPACK_ROOT}"
    "lib"
    )

  if(NOT _la_lapacke_bundled)
    # well, it's not bundled, let's give it another shot by also linking the
    # generic liblapacke.so
    check_cblas_lapacke(_la_lapacke_library
      LAPACKE
      LAPACKE_cheev
      ${_include_dirs}
      "lapacke.h"
      ""
      "lapacke;${_lapack_libs}"
      ""
      "${BLAS_ROOT};${LAPACK_ROOT}"
      "lib"
      )

    if(_la_lapacke_library)
      list(PREPEND LAPACK_LIBRARIES ${_la_lapacke_library})
      unset(_la_lapacke_library)
      unset(_include_dirs)
    else()
      message(FATAL_ERROR "Cannot find a working LAPACKE library to link to!")
    endif()
  endif()
endif()

set(CMAKE_REQUIRED_QUIET OFF)

# create target
# adapted from: https://github.com/Kitware/CMake/blob/master/Modules/FindLAPACK.cmake#L588-L603
if(BLAS_FOUND AND LAPACK_FOUND AND NOT TARGET Math::LA)
  message(STATUS "Using ${VLX_LA_VENDOR} as linear algebra backend")

  add_library(Math::LA INTERFACE IMPORTED)

  target_compile_definitions(Math::LA
    INTERFACE
      "${_la_compile_definitions}"
    )

  if(_la_compiler_flags)
    target_compile_options(Math::LA
      INTERFACE
        "${_la_compiler_flags}"
      )
  endif()

  if(_la_include_dirs)
    # not strictly necessary, but clean up duplicates
    list(REMOVE_DUPLICATES _la_include_dirs)

    message(STATUS "BLAS/LAPACK include directories: ${_la_include_dirs}")

    target_include_directories(Math::LA
      SYSTEM
      INTERFACE
        "${_la_include_dirs}"
      )
  endif()

  # not strictly necessary, but clean up duplicates
  list(REMOVE_DUPLICATES LAPACK_LIBRARIES)

  message(STATUS "BLAS/LAPACK libraries: ${LAPACK_LIBRARIES}")

  target_link_libraries(Math::LA
    INTERFACE
      "${LAPACK_LIBRARIES}"
  )

  unset(_la_compile_definitions)
  unset(_la_compiler_flags)
  unset(_la_include_dirs)
  unset(_la_libraries)
endif()
