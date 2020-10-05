#.rst:
#
# Enables architecture-specific compiler flags.
#
# Variables used::
#
#   ENABLE_ARCH_FLAGS
#
# autocmake.yml configuration::
#
#   docopt: "--arch-flags=<ARCH_FLAGS> Enable architecture-specific compiler flags [default: True]."
#   define: "'-DENABLE_ARCH_FLAGS={0}'.format(arguments['--arch-flags'])"

option(ENABLE_ARCH_FLAGS "Enable architecture-specific compiler flags" ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

if(ENABLE_ARCH_FLAGS)
  include(${CMAKE_CURRENT_LIST_DIR}/set_compiler_flag.cmake)
  # iterate over list of flags and use the first one that is compatible with the
  # compiler in use
  set_compiler_flag(_arch_flag
    FLAGS
      "-march=native"  # valid wth GNU and probably Clang too
      "-xHost"  # valid with Intel
    )
  message(STATUS "Adding architecture-specific compiler flag: ${_arch_flag}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_arch_flag}")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/GNU.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Intel.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Clang.CXX.cmake)

string(REPLACE " " ";" _cmake_cxx_flags ${CMAKE_CXX_FLAGS})
string(REPLACE " " ";" _cmake_cxx_flags_release ${CMAKE_CXX_FLAGS_RELEASE})
string(REPLACE " " ";" _cmake_cxx_flags_relwithdebinfo ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
