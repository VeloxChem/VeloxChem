# Adapted from
# https://github.com/dev-cafe/cmake-cookbook/blob/master/chapter-07/recipe-03/c-cxx-example/set_compiler_flag.cmake

# Given a list of flags, this stateless function will try each, one at a time,
# and set result to the first flag that works.
# If none of the flags works, result is "".
# If the REQUIRED key is given and no flag is found, a FATAL_ERROR is raised.
#
# Call is:
# set_compiler_flag(result <REQUIRED> FLAGS flag1 flag2 ...)
#
# Example:
# set_compiler_flag(working_compile_flag REQUIRED FLAGS "-Wall" "-warn all")

include(CheckCXXCompilerFlag)

function(set_compiler_flag _result)
  set(options REQUIRED)
  set(oneValueArgs)
  set(multiValueArgs FLAGS)
  cmake_parse_arguments(_flagger
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  # Silently check compiler flags
  set(restore_CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET})
  set(CMAKE_REQUIRED_QUIET TRUE)

  # loop over all flags, try to find the first which works
  set(_flag_found FALSE)
  foreach(flag IN ITEMS ${_flagger_FLAGS})
    unset(_flag_works CACHE)
    check_cxx_compiler_flag("${flag}" _flag_works)

    # if the flag works, use it, and exit
    # otherwise try next flag
    if(_flag_works)
      set(${_result} "${flag}" PARENT_SCOPE)
      set(_flag_found TRUE)
      break()
    endif()
  endforeach()

  # raise an error if no flag was found
  if(_flagger_REQUIRED AND NOT _flag_found)
    message(FATAL_ERROR "None of the required flags were supported")
  endif()

  # Restore CMAKE_REQUIRED_QUIET
  set(CMAKE_REQUIRED_QUIET ${restore_CMAKE_REQUIRED_QUIET})
endfunction()
