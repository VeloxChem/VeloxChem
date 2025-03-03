target_sources(vlxobjs
  PRIVATE
    RIFockDriver.cpp
    T3CDistributor.cpp
    T3CUtils.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
