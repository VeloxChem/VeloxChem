target_sources(vlxobjs
  PRIVATE
    T2CGeom10SumDistributor.cpp
    T2CGeom10SumTwoDistributor.cpp
    T2CUtils.cpp
    DenseMatrixDistributor.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
