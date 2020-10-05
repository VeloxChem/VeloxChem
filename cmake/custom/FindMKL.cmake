# we mimick the include directories, compiler options, and link line suggested by the advisor:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
add_library(MKL INTERFACE)

# locate mkl.h
find_path(MKL_INCLUDE_DIRS
  NAMES
    mkl.h
  HINTS
    ${CMAKE_INSTALL_PREFIX}/include
  )

# locate MKL libraries
# With this list we choose to **always** link against the Intel OpenMP runtime,
# rather than the one provided by GNU
set(MKL_LIBRARIES)
foreach(_l IN ITEMS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
  find_library(_x
    NAMES
      ${_l}
    HINTS
      ${CMAKE_INSTALL_PREFIX}/lib
    )
  list(APPEND MKL_LIBRARIES ${_x})
  unset(_x CACHE)
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL
  DEFAULT_MSG
  MKL_INCLUDE_DIRS
  MKL_LIBRARIES
  )

# <<< create a proper CMake target >>>
# compile options
target_compile_options(MKL
  INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:AppleClang>>:-m64>
  )

# include directories
target_include_directories(MKL
  INTERFACE
    ${MKL_INCLUDE_DIRS}
  )

# link libraries
find_package(Threads QUIET)
target_link_libraries(MKL
  INTERFACE
    ${MKL_LIBRARIES}
    Threads::Threads
    m
    dl
  )

# TODO One needs to call mkl_set_interface_layer and mkl_set_threading_layer
# to set up runtime properly. Since these things are known when configuring
# the project, it would be nice to generate a source file making the calls.
