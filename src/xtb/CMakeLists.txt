find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
  pkg_search_module(xTB
    IMPORTED_TARGET
    xtb
    )
  
  if(TARGET PkgConfig::xTB) 
    message(STATUS "Found xTB ${xTB_VERSION}: ${xTB_LINK_LIBRARIES}")
    set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/XTBDriver.cpp
      PROPERTIES
        COMPILE_DEFINITIONS ENABLE_XTB
      )
    target_link_libraries(vlxobjs
      PUBLIC
        PkgConfig::xTB
      )
  endif()
else()
   message(STATUS "Missing pkg-config: cannot attempt finding xTB!")
endif()

target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/XTBDriver.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
