target_sources(vlxobjs
  PRIVATE
    LebedevLaikovQuadrature.cpp
    Log3Quadrature.cpp
    DensityGridGenerator.cpp
    DftSubMatrix.cpp
    GridBox.cpp
    GridDriver.cpp
    GridPartitioner.cpp
    MolecularGrid.cpp
    PartitionFunc.cpp
    XCIntegrator.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
