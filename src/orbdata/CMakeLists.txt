target_sources(vlxobjs
  PRIVATE
    AODensityMatrix.cpp
    AtomBasis.cpp
    GtoFunc.cpp
    MolecularBasis.cpp
    GtoBlock.cpp
    BasisFunction.cpp
    GtoPairBlock.cpp 	
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
