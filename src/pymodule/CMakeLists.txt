if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  set(_mkl_interface_layer_ "GNU,LP64")
  set(_mkl_threading_layer_ "GNU")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(_mkl_interface_layer_ "LP64")
  set(_mkl_threading_layer_ "INTEL")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(_mkl_interface_layer_ "LP64")
  set(_mkl_threading_layer_ "INTEL")
endif()

# generate MKL runtime configuration function
configure_file(
  mklconf.py.in
  ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/mklconf.py
  @ONLY
  )

# we glob the Python files in src/pymodule and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
  GLOB
    _veloxchemlib_pys
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py
  )

# link the Python files under the build folder
foreach(_py IN LISTS _veloxchemlib_pys)
  get_filename_component(__py ${_py} NAME)
  file(
    CREATE_LINK
      ${_py}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${__py}
    COPY_ON_ERROR
    SYMBOLIC
    )
endforeach()

list(APPEND _veloxchemlib_pys
  ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/mklconf.py
  )

install(
  FILES
    ${_veloxchemlib_pys}
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
