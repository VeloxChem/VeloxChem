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
# note that __init__.py is skipped
foreach(_py_path IN LISTS _veloxchemlib_pys)
  get_filename_component(_py_name ${_py_path} NAME)
  if(_py_name STREQUAL "__init__.py")
    continue()
  endif()
  file(
    CREATE_LINK
      ${_py_path}
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/${_py_name}
    COPY_ON_ERROR
    SYMBOLIC
    )
endforeach()

# process commit hash
set(_vlx_commit "unknown")

if(NOT DEFINED ENABLE_COMMIT_HASH)
  set(ENABLE_COMMIT_HASH FALSE)
endif()

if(ENABLE_COMMIT_HASH)
  find_program(_git_exe git)
  if(_git_exe)
    execute_process(
      COMMAND "${_git_exe}" rev-parse --short HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE _commit_hash_str
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE _git_ret
    )
    if(_git_ret EQUAL 0 AND NOT _commit_hash_str STREQUAL "")
      set(_vlx_commit "${_commit_hash_str}")
      message(STATUS "Setting _vlx_commit: ${_commit_hash_str}")
      file(REMOVE
          ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py)
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
          DESTINATION ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR})
      file(APPEND
          ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
          "\n__commit__ = \"${_vlx_commit}\"\n")
      list(APPEND
        _veloxchemlib_pys
        ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py)          
    endif()
  endif()
endif()

if(${_vlx_commit} STREQUAL "unknown")
  file(REMOVE
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py)
  file(CREATE_LINK
      ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
      ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
      COPY_ON_ERROR
      SYMBOLIC)
  list(APPEND
      _veloxchemlib_pys
      ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py)             
endif()

install(
  FILES
    ${_veloxchemlib_pys}
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
