target_sources(vlxobjs
  PRIVATE
    LebedevLaikovQuadrature.cpp
    Log3Quadrature.cpp
    GridBox.cpp
    GridDriver.cpp
    GridPartitioner.cpp
    MolecularGrid.cpp
    PartitionFunc.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
