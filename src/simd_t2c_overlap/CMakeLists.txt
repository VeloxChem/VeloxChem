target_sources(vlxobjs
  PRIVATE
    SimdOverlapDriver.cpp
    SimdOverlapFunc.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
