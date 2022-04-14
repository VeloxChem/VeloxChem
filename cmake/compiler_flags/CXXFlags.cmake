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

# code needs C++17 at least
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# do not use compiler extensions to the C++ standard
set(CMAKE_CXX_EXTENSIONS FALSE)
# generate a JSON database of compiler commands (useful for LSP IDEs)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

if(ENABLE_ARCH_FLAGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(_arch_flag "-march=native")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if(WIN32) # use AVX2 on Windows
      set(_arch_flag "/arch:AVX2")
    else()
      set(_arch_flag "-march=native")
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(_arch_flag "-xHost")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_arch_flag}")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/GNU.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Intel.CXX.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Clang.CXX.cmake)

string(REPLACE " " ";" _cmake_cxx_flags ${CMAKE_CXX_FLAGS})
string(REPLACE " " ";" _cmake_cxx_flags_release ${CMAKE_CXX_FLAGS_RELEASE})
string(REPLACE " " ";" _cmake_cxx_flags_relwithdebinfo ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
