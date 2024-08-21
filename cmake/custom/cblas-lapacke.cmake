function(CHECK_CBLAS_LAPACKE 
  LIBRARIES
  _prefix
  _name
  _includedir
  _header
  _flags 
  _list 
  _deps 
  _addlibdir 
  _subdirs
  )
  # This function checks for the existence of the combination of libraries
  # given by _list.  If the combination is found, this checks whether can link
  # against that library combination using the name of a routine given by _name
  # including the header given by _header in _includedir and using the linker flags given by _flags.  If the combination of libraries is
  # found and passes the link test, ${LIBRARIES} is set to the list of complete
  # library paths that have been found.  Otherwise, ${LIBRARIES} is set to FALSE.

  set(_libraries_work TRUE)
  set(_libraries)
  set(_combined_name)

  if(BLA_STATIC)
    if(WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else()
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif()
  else()
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # for ubuntu's libblas3gf and liblapack3gf packages
      set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
    endif()
  endif()

  set(_extaddlibdir "${_addlibdir}")
  if(WIN32)
    list(APPEND _extaddlibdir ENV LIB)
  elseif(APPLE)
    list(APPEND _extaddlibdir ENV DYLD_LIBRARY_PATH)
  else()
    list(APPEND _extaddlibdir ENV LD_LIBRARY_PATH)
  endif()
  list(APPEND _extaddlibdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")

  foreach(_library ${_list})
    if(_library MATCHES "^-")
      # Respect linker flags as-is (required by MKL)
      list(APPEND _libraries "${_library}")
    else()
      string(REGEX REPLACE "[^A-Za-z0-9]" "_" _lib_var "${_library}")
      string(APPEND _combined_name "_${_lib_var}")
      if(NOT "${_deps}" STREQUAL "")
        string(APPEND _combined_name "_deps")
      endif()
      if(_libraries_work)
        find_library(${_prefix}_${_lib_var}_LIBRARY
          NAMES ${_library}
          NAMES_PER_DIR
          PATHS ${_extaddlibdir}
          PATH_SUFFIXES ${_subdirs}
        )
        mark_as_advanced(${_prefix}_${_lib_var}_LIBRARY)
        list(APPEND _libraries ${${_prefix}_${_lib_var}_LIBRARY})
        set(_libraries_work ${${_prefix}_${_lib_var}_LIBRARY})
      endif()
    endif()
  endforeach()

  foreach(_flag ${_flags})
    string(REGEX REPLACE "[^A-Za-z0-9]" "_" _flag_var "${_flag}")
    string(APPEND _combined_name "_${_flag_var}")
  endforeach()
  if(_libraries_work)
    # Test this combination of libraries.
    set(CMAKE_REQUIRED_INCLUDES ${_includedir})
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${_libraries} ${_deps})
    include(CheckSymbolExists)
    check_symbol_exists("${_name}" "${_header}" ${_prefix}${_combined_name}_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif()

  if(_libraries_work)
    if("${_list}" STREQUAL "")
      set(_libraries "${LIBRARIES}-PLACEHOLDER-FOR-EMPTY-LIBRARIES")
    else()
      list(APPEND _libraries ${_deps})
    endif()
  else()
    set(_libraries FALSE)
  endif()
  set(${LIBRARIES} "${_libraries}" PARENT_SCOPE)
endfunction()
