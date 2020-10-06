# Export compile commands for each file to JSON
# This is useful for static analysis tools and linters
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED COMPONENTS CXX)

find_package(OpenMP REQUIRED COMPONENTS CXX)

include(${CMAKE_CURRENT_LIST_DIR}/FindMKL.cmake)

# figure out where to put the Python module
if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within CMAKE_INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
# install Python module under CMAKE_INSTALL_LIBDIR
# if that is "lib64", the use just "lib"
set(_lib "${CMAKE_INSTALL_LIBDIR}")
if(CMAKE_INSTALL_LIBDIR STREQUAL "lib64")
  set(_lib "lib")
endif()
file(TO_NATIVE_PATH "${_lib}/${PYMOD_INSTALL_LIBDIR}/veloxchem" PYMOD_INSTALL_FULLDIR)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})

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
  RELATIVE
    ${PROJECT_SOURCE_DIR}/src/pymodule
  ${PROJECT_SOURCE_DIR}/src/pymodule/*.py
  )

# link the Python files under the build folder
foreach(_f IN LISTS _veloxchemlib_pys)
  file(
    CREATE_LINK
      ${PROJECT_SOURCE_DIR}/src/pymodule/${_f}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${_f}
    COPY_ON_ERROR
    SYMBOLIC
    )
endforeach()

add_subdirectory(src)

# handle folder with basis sets
# 1. symlink under the build tree
file(
  CREATE_LINK
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/basis
  COPY_ON_ERROR
  SYMBOLIC
  )
# 2. install rules for basis sets folder
install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )

enable_testing()
include(CTest)
# This must come last!!
add_subdirectory(unit_tests)
add_subdirectory(python_tests)
