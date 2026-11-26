#.rst:
#
# The hardware specific math library of the dense linear algebra: Accelerate
# on Apple systems, and MKL or OpenBLAS elsewhere through FindLAPACK, which
# finds BLAS as well. When disabled, Eigen is used. The routines are the
# Fortran symbols of the LP64 interface.
#
# Variables used::
#
#   VLX_USE_MATHLIB
#   VLX_MATH_LIBRARY
#   BLA_VENDOR
#
# Imported targets::
#
#   Math::Mathlib
#     The math library to link, always defined.

option_with_print(VLX_USE_MATHLIB "Link a hardware specific math library for dense linear algebra" OFF)

# always defined, so that the sources link it unconditionally
add_library(Math::Mathlib INTERFACE IMPORTED)

if(VLX_USE_MATHLIB)

  target_compile_definitions(Math::Mathlib
    INTERFACE
      VLX_USE_MATHLIB
    )

  if(APPLE)

    target_compile_definitions(Math::Mathlib
      INTERFACE
        ACCELERATE_NEW_LAPACK
      )

    target_link_libraries(Math::Mathlib
      INTERFACE
        "-framework Accelerate"
      )

    message(STATUS "Using the Accelerate framework as the math library")

  elseif(DEFINED VLX_MATH_LIBRARY AND NOT "${VLX_MATH_LIBRARY}" STREQUAL "")

    message(STATUS "Using the math library set by VLX_MATH_LIBRARY: ${VLX_MATH_LIBRARY}")

    target_link_libraries(Math::Mathlib
      INTERFACE
        ${VLX_MATH_LIBRARY}
      )

  else()

    find_package(LAPACK REQUIRED)

    target_link_libraries(Math::Mathlib
      INTERFACE
        LAPACK::LAPACK
      )

    get_property(_mathlib_libraries TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES)
    message(STATUS "Using the LAPACK found by CMake: ${_mathlib_libraries}")

  endif()

endif()
