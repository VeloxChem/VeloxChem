target_sources(vlxobjs
  PRIVATE
    SimdThreeCenterElectronRepulsionDriver.cpp
    SimdThreeCenterElectronRepulsionFunc.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
