target_sources(vlxobjs
  PRIVATE
    RIFockDriver.cpp
    RIFockGradDriver.cpp
    T3CDistributor.cpp
    T3CGeom010SumDistributor.cpp
    T3CGeom0X0Distributor.cpp
    T3CGeom100SumDistributor.cpp
    T3CGeomX00Distributor.cpp
    T3CUtils.cpp
    T3CLocalDistributor.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
