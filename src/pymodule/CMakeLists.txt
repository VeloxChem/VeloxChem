# Install Python files from the source tree directly.
# install(DIRECTORY ...) preserves subdirectory structure (e.g. mofbuilder/)
# and resolves file matching at install time, so new or deleted .py files are
# automatically picked up without a CMake reconfigure.
# __init__.py is excluded from the directory install to allow commit-hash
# injection and is handled separately below.
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${PYMOD_INSTALL_FULLDIR}
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "__init__.py" EXCLUDE
  PATTERN "__pycache__" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "*.egg-info" EXCLUDE
  )

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
    endif()
  endif()
endif()

if(NOT _vlx_commit STREQUAL "unknown")
  message(STATUS "Setting __commit__: ${_vlx_commit}")
  file(READ ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py _init_content)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/veloxchem_init.py
       "${_init_content}\n__commit__ = \"${_vlx_commit}\"\n")
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/veloxchem_init.py
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
    RENAME __init__.py
  )
else()
  install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
    DESTINATION ${PYMOD_INSTALL_FULLDIR}
  )
endif()
