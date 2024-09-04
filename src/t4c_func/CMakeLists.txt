target_sources(vlxobjs
  PRIVATE
    FockDriver.cpp
    T4CLocalDistributor.cpp
    T4CMatrixDistributor.cpp
    T4CUtils.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
