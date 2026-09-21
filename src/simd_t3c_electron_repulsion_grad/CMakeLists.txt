target_sources(vlxobjs
  PRIVATE
    SimdThreeCenterElectronRepulsionGradientDriver.cpp
    SimdThreeCenterElectronRepulsionGradientRows010.cpp
    SimdThreeCenterElectronRepulsionGradientRows100.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
