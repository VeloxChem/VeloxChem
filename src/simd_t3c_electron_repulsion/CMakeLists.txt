target_sources(vlxobjs
  PRIVATE
    SimdThreeCenterElectronRepulsionCtrVrrSSS.cpp
    SimdThreeCenterElectronRepulsionDriver.cpp
    SimdThreeCenterElectronRepulsionFunc.cpp
    SimdThreeCenterElectronRepulsionRecSSS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
