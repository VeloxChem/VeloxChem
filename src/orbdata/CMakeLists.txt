target_sources(vlxobjs
  PRIVATE
    AtomBasis.cpp
    BasisFunction.cpp
    GtoBlock.cpp
    GtoFunc.cpp
    MolecularBasis.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
