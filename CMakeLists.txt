# set minimum cmake version
cmake_minimum_required(VERSION 3.16)

# project name
project(VeloxChem LANGUAGES CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# if CMAKE_BUILD_TYPE undefined, we set it to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Options handling utilities
include(CMakeDependentOption)
# Macro for printing an option in a consistent manner
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# Syntax: print_option(<option to print> <was specified>)
macro(print_option variable default)
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
  else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
  endif()
endmacro()

# Wraps an option with default ON/OFF. Adds nice messaging to option()
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# Syntax: option_with_print(<option name> <description> <default value>)
macro(option_with_print variable msge default)
  print_option(${variable} ${default})
  option(${variable} ${msge} ${default})
endmacro()

# Wraps an option with a default other than ON/OFF and prints it
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard)
# NOTE
#  1. Can't combine with above b/c CMake handles ON/OFF options specially
#  2. CMake variables are always defined so need to further check for if
#     they are the NULL string. This is also why we need the force
# Syntax: option_with_default(<option name> <description> <default value>)
macro(option_with_default variable msge default)
  print_option(${variable} "${default}")
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    set(${variable} "${default}" CACHE STRING ${msge} FORCE)
  endif()
endmacro()

# included cmake modules
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_cxx.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/compiler_flags/CXXFlags.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_default_build_paths.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_safeguards.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/custom/alignment.cmake)

# Currently, scikit-build does not support FindPython, so we convert the
# provided hints ourselves.
if(SKBUILD)
  set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
  set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
  set(Python_LIBRARY "${PYTHON_LIBRARY}")
  set(DUMMY "${PYTHON_VERSION_STRING}")  # Not needed, silences a warning
endif()
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter Development NumPy)
# download pybind11 if not available on the system
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_pybind11.cmake)
# TODO code coverage
include(${PROJECT_SOURCE_DIR}/cmake/custom/vlx.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/downloaded/autocmake_save_flags.cmake)
