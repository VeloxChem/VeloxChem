target_sources(vlxobjs
  PRIVATE
    PairDensityBecke88.cpp
    PairDensityGridGenerator.cpp
    PairDensityLYP.cpp
    PairDensityPBE_C.cpp
    PairDensityPBE_X.cpp
    PairDensitySlater.cpp
    PairDensityVWN.cpp
    XCIntegratorForPDFT.cpp
    XCPairDensityFunctional.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
