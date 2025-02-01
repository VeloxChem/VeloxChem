# TODO: move dftd4 search function to cmake/custom/

set(_dftd4_inc_search_paths)
if(DEFINED ENV{DFTD4_HOME} OR DEFINED ENV{DFTD4_INCLUDE_DIR})
    list(APPEND _dftd4_inc_search_paths PATHS)
    if(DEFINED ENV{DFTD4_HOME})
        list(APPEND _dftd4_inc_search_paths $ENV{DFTD4_HOME}/include)
    endif()
    if(DEFINED ENV{DFTD4_INCLUDE_DIR})
        list(APPEND _dftd4_inc_search_paths $ENV{DFTD4_INCLUDE_DIR})
    endif()
endif()

# Note: only search for dftd4 when DFTD4_HOME or DFTD4_INCLUDE_DIR is defined
if(_dftd4_inc_search_paths)
    message(STATUS "Searching for DFTD4...")
    find_path(_dftd4_include_dirs
        NAMES dftd4.h
        ${_dftd4_inc_search_paths}
        HINTS ${CMAKE_PREFIX_PATH}/include
    )
endif()

set(_dftd4_lib_search_paths)
if(DEFINED ENV{DFTD4_HOME} OR DEFINED ENV{DFTD4_LIB_DIR})
    list(APPEND _dftd4_lib_search_paths PATHS)
    if(DEFINED ENV{DFTD4_HOME})
        list(APPEND _dftd4_lib_search_paths $ENV{DFTD4_HOME}/lib)
    endif()
    if(DEFINED ENV{DFTD4_LIB_DIR})
        list(APPEND _dftd4_lib_search_paths $ENV{DFTD4_LIB_DIR})
    endif()
endif()

# Note: only search for dftd4 when DFTD4_HOME or DFTD4_LIB_DIR is defined
if(_dftd4_lib_search_paths)
    find_library(_dftd4_libraries
        NAMES "dftd4"
        ${_dftd4_lib_search_paths}
        HINTS ${CMAKE_PREFIX_PATH}/lib
    )
endif()

if(_dftd4_inc_search_paths AND _dftd4_lib_search_paths)
  if(_dftd4_include_dirs AND _dftd4_libraries)
    message(STATUS "Found DFTD4: ${_dftd4_libraries}")
    target_compile_definitions(vlxobjs
      PUBLIC
        ENABLE_DFTD4
      )
    target_include_directories(vlxobjs
      SYSTEM
      PRIVATE
        ${_dftd4_include_dirs}
      )
    target_link_libraries(vlxobjs
      PUBLIC
        ${_dftd4_libraries}
      )
  else()
    message(STATUS "Could not find DFTD4.")
  endif()
endif()

target_sources(vlxobjs
  PRIVATE
    AtomicPartialChargesModel.cpp
    AtomicRadii.cpp
    ChemicalElement.cpp
    DispersionModel.cpp
    Molecule.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
