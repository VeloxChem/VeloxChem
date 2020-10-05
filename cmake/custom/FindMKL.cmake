# we mimick the include directories, compiler options, and link line suggested by the advisor:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
add_library(MKL INTERFACE)

# locate mkl.h
find_path(_mkl_h
  NAMES
    mkl.h
  HINTS
    ${CMAKE_INSTALL_PREFIX}/include
  )

# locate MKL single dynamic library (mkl_rt)
find_library(_mkl_libs
  NAMES
    mkl_rt
  HINTS
    ${CMAKE_INSTALL_PREFIX}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG _mkl_h _mkl_libs)

# <<< create a proper CMake target >>>

# compile options
target_compile_options(MKL
  INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:AppleClang>>:-m64>
  )

# include directories
target_include_directories(MKL
  INTERFACE
    ${_mkl_h}
  )

# link libraries
find_package(Threads QUIET)
target_link_libraries(MKL
  INTERFACE
    ${_mkl_libs}
    Threads::Threads
    m
    dl
  )

# TODO One needs to call mkl_set_interface_layer and mkl_set_threading_layer
# to set up runtime properly. Since these things are known when configuring
# the project, it would be nice to generate a source file making the calls.
