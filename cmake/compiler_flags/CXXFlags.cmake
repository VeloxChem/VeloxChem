#.rst:
#
# Manages C++ compiler flags.
#
# There is one user-facing option to enable architecture-specific compiler
# flags.
# The complete list of flags is built as:
#
#   CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_<CONFIG> ARCH_FLAG VLX_CXX_FLAGS EXTRA_CXXFLAGS
#
# where:
#
# - ``CMAKE_CXX_FLAGS`` is initialized by the contents of the ``CXXFLAGS``
#   environment variable when configuring. Default is empty.
# - ``CMAKE_CXX_FLAGS_<CONFIG>`` are build-type specific compiler flags.
#   The defaults are compiler-dependent: have a look at the ``GNU.CXX.cmake``,
#   ``Clang.CXX.cmake``, and ``Intel.CXX.cmake`` files.
# - ``ARCH_FLAG`` is the architecture-dependent optimization flag, *e.g.*
#   vectorization. Default is empty.
# - ``VECLIB_FLAG`` selects the vector math library, so that the vectorizer can
#   replace the scalar math calls of the integral kernels, the exponential above
#   all, with their vector counterparts. Default is empty.
# - ``VLX_CXX_FLAGS`` are VeloxChem-specific flags to be used for all builds.
#   The defaults are compiler-dependent: have a look at the ``GNU.CXX.cmake``,
#   ``Clang.CXX.cmake``, and ``Intel.CXX.cmake`` files.
# - ``EXTRA_CXXFLAGS`` useful if you need to append certain flags to the full
#   list, *e.g.* to override previous compiler flags without touching the CMake
#   scripts.  Default is empty.
#
# Variables used::
#
#   ENABLE_ARCH_FLAGS
#   ENABLE_VECLIB
#   EXTRA_CXXFLAGS
#
# Variables modified::
#
#   CMAKE_CXX_FLAGS
#
# Environment variables used::
#
#   CXXFLAGS

option_with_print(ENABLE_ARCH_FLAGS "Enable architecture-specific compiler flags" ON)
option_with_print(ENABLE_VECLIB "Enable the vector math library of the platform" ON)

# code needs C++20 at least
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# do not use compiler extensions to the C++ standard
set(CMAKE_CXX_EXTENSIONS FALSE)
# generate a JSON database of compiler commands (useful for LSP IDEs)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
# position-independent code
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
# visibility levels
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)

set(ARCH_FLAG "")
if(ENABLE_ARCH_FLAGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(ARCH_FLAG "-march=native")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if(WIN32) # use AVX2 on Windows
      set(ARCH_FLAG "/arch:AVX2")
    else()
      set(ARCH_FLAG "-march=native")
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(ARCH_FLAG "-xHost")
  endif()
endif()

# The vector math library is selected by the platform: Apple systems carry the
# two wide double precision routines in libsystem_m, while glibc carries them in
# libmvec. The probe compiles and links a vectorized loop, as a flag which the
# driver accepts still leaves undefined symbols when the library is absent.

set(VECLIB_FLAG "")
if(ENABLE_VECLIB AND CMAKE_CXX_COMPILER_ID MATCHES Clang AND NOT WIN32)
  if(APPLE)
    set(_veclib "Darwin_libsystem_m")
  else()
    set(_veclib "libmvec")
  endif()

  include(CheckCXXSourceCompiles)

  set(_saved_flags "${CMAKE_REQUIRED_FLAGS}")
  set(_saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -O2 -fveclib=${_veclib} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} OpenMP::OpenMP_CXX)

  check_cxx_source_compiles("
    #include <cmath>
    int main()
    {
        double x[64], y[64];
        for (int i = 0; i < 64; i++) x[i] = -0.5 * i;
    #pragma omp simd
        for (int i = 0; i < 64; i++) y[i] = std::exp(x[i]);
        return (y[0] > 0.5) ? 0 : 1;
    }" VECLIB_LINKS)

  set(CMAKE_REQUIRED_FLAGS "${_saved_flags}")
  set(CMAKE_REQUIRED_LIBRARIES "${_saved_libs}")

  if(VECLIB_LINKS)
    set(VECLIB_FLAG "-fveclib=${_veclib}")
  endif()
endif()
message(STATUS "Vector math library flag: ${VECLIB_FLAG}")

set(VLX_CXX_FLAGS "")
include(${CMAKE_CURRENT_LIST_DIR}/GNU.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Intel.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Clang.CXX.cmake)
