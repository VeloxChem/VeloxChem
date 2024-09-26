target_sources(vlxobjs
  PRIVATE
    CubicGrid.cpp
    VisualizationDriver.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
