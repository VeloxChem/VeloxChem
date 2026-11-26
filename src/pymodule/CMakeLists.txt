# Install Python files from the source tree directly.
# install(DIRECTORY ...) preserves subdirectory structure (e.g. mofbuilder/)
# and resolves file matching at install time, so new or deleted .py files are
# automatically picked up without a CMake reconfigure.
# All __init__.py files are installed here (e.g. mofbuilder/__init__.py). The
# top-level __init__.py is replaced below when a commit hash is injected.
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${PYMOD_INSTALL_FULLDIR}
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "__pycache__" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "*.egg-info" EXCLUDE
  )

# Optionally append the commit hash to the top-level __init__.py.
# The hash is resolved at install time (not configure time) so that
# re-installing from an existing build directory picks up new commits
# without a CMake reconfigure.
if(NOT DEFINED ENABLE_COMMIT_HASH)
  set(ENABLE_COMMIT_HASH FALSE)
endif()

if(ENABLE_COMMIT_HASH)
  find_program(_git_exe git)
  if(_git_exe)
    install(CODE "
      execute_process(
        COMMAND \"${_git_exe}\" rev-parse --short HEAD
        WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\"
        OUTPUT_VARIABLE _vlx_commit
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _git_ret
      )
      if(_git_ret EQUAL 0 AND NOT _vlx_commit STREQUAL \"\")
        file(READ \"${CMAKE_CURRENT_SOURCE_DIR}/__init__.py\" _init_content)
        file(WRITE \"${CMAKE_CURRENT_BINARY_DIR}/veloxchem_init.py\"
             \"\${_init_content}\\n__commit__ = \\\"\${_vlx_commit}\\\"\\n\")
        # Remove the plain file installed above so that the commit hash is
        # always injected, regardless of file timestamps.
        file(REMOVE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR}/__init__.py\")
        file(INSTALL
             DESTINATION \"\${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR}\"
             TYPE FILE RENAME \"__init__.py\"
             FILES \"${CMAKE_CURRENT_BINARY_DIR}/veloxchem_init.py\")
        message(STATUS \"Setting __commit__: \${_vlx_commit}\")
      endif()
    ")
  endif()
endif()
