target_sources(vlxobjs
  PRIVATE
    AtomBasis.cpp
    GtoFunc.cpp
    MolecularBasis.cpp
    GtoBlock.cpp
    BasisFunction.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
