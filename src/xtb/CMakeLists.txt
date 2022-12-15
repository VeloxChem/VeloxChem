if(DEFINED ENV{XTBHOME})
  message(STATUS "Searching for xTB in $ENV{XTBHOME}")

  find_path(_xtb_include_dirs
    NAMES
      xtb.h
    PATHS ENV
      XTBHOME
    PATH_SUFFIXES
      include
    )

  find_library(_xtb_libraries
    NAMES
      "xtb"
    PATHS ENV
      XTBHOME
    PATH_SUFFIXES
      lib
    )

  try_run(_xtb_RUNS _xtb_COMPILES
    ${CMAKE_CURRENT_BINARY_DIR}/xtb_discovery
    SOURCES
      "${PROJECT_SOURCE_DIR}/cmake/custom/xtb_discovery.cpp"
    LINK_LIBRARIES
      "${_xtb_libraries}"
    ARGS
      "-DINCLUDE_DIRECTORIES:STRING=${_xtb_include_dirs}"
    )

  if(_xtb_COMPILES AND (_xtb_RUNS EQUAL 0))
    message(STATUS "Found xTB: ${_xtb_libraries}")

    target_compile_definitions(vlxobjs
      PUBLIC
        ENABLE_XTB
      )

    target_include_directories(vlxobjs
      SYSTEM
      PRIVATE
        ${_xtb_include_dirs}
      )

    target_link_libraries(vlxobjs
      PUBLIC
        ${_xtb_libraries}
      )
  else()
    message(STATUS "Could not find a working xTB installation!")
  endif()
endif()

target_sources(vlxobjs
  PRIVATE
    XtbDriver.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
