target_sources(vlxobjs
  PRIVATE
    ThreeCenterOverlapDriver.cpp
    ThreeCenterOverlapPrimRecSS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
