target_sources(vlxobjs
  PRIVATE
    SimdBoysFunc.cpp
    SimdBoysTable.cpp
    SimdCoordinates.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
