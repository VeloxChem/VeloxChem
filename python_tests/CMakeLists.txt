# we glob the Python test files and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
  GLOB
    _vlx_pytests
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  RELATIVE
    ${CMAKE_CURRENT_LIST_DIR}
  test_*.py
  )
# append conftest.py to the list
list(APPEND _vlx_pytests conftest.py)

foreach(_f IN LISTS _vlx_pytests)
  # link the Python files under the build folder
  file(
    CREATE_LINK
      ${CMAKE_CURRENT_LIST_DIR}/${_f}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../tests/${_f}
    COPY_ON_ERROR
    SYMBOLIC
    )

  if(${_f} MATCHES "conftest.py")
    break()
  endif()

  # extract the name and label for the test from the file name
  # by removing the "test_" prefix and the ".py" suffix
  string(REGEX REPLACE "test_(.*).py" "\\1" _name ${_f})
  add_test(
    NAME
      ${_name}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python_EXECUTABLE} -m pytest -v -s ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../tests/${_f}
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/..
    )

  set_tests_properties(${_name}
    PROPERTIES
      LABELS
        "integration;${_name}"
    )
endforeach()

# handle folder with reference input files
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../tests/inputs)
# we glob the reference input files in inputs and let CMake add a rule such that
# the glob is repeated every time we rebuild.
# This is NOT RECOMMENDED by CMake
# (https://cmake.org/cmake/help/v3.16/command/file.html#filesystem) but you only
# live once!
file(
  GLOB
    _ref_inputs
  LIST_DIRECTORIES
    FALSE
  CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_LIST_DIR}/inputs/*
  )
# 1. symlink under the build tree
foreach(_input IN LISTS _ref_inputs)
 get_filename_component(__input ${_input} NAME)
 file(
   CREATE_LINK
     ${_input}
     ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/../tests/inputs/${__input}
   COPY_ON_ERROR
   SYMBOLIC
   )
endforeach()

# install the current folder
install(
  DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}/
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}/tests
  PATTERN "py_plot_comb.py" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "Makefile" EXCLUDE
  )
