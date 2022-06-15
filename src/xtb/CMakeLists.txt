find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
  pkg_search_module(xTB
    IMPORTED_TARGET
    xtb
    )
  
  if(TARGET PkgConfig::xTB) 
    message(STATUS "Found xTB ${xTB_VERSION}: ${xTB_LINK_LIBRARIES}")
    target_compile_definitions(vlxobjs
      PUBLIC
        ENABLE_XTB
      )
    target_link_libraries(vlxobjs
      PUBLIC
        PkgConfig::xTB
      )
  else()
    message(STATUS "Could not find xTB!")
  endif()
else()
   message(STATUS "Missing pkg-config: cannot attempt finding xTB!")
endif()

target_sources(vlxobjs
  PRIVATE
    XTBDriver.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
