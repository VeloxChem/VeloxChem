target_sources(vlxobjs
  PRIVATE
    SimdTwoCenterElectronRepulsionDriver.cpp
    SimdTwoCenterElectronRepulsionFunc.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
