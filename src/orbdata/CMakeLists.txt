target_sources(vlxobjs
  PRIVATE
    AODensityMatrix.cpp
    AOIndices.cpp
    AtomBasis.cpp
    BasisFunction.cpp
    BlockedGtoPairBlock.cpp
    GtoBlock.cpp
    GtoFunc.cpp
    GtoPairBlock.cpp
    GtoPairBlockFunc.cpp
    MolecularBasis.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
