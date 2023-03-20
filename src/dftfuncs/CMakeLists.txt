target_sources(vlxobjs
  PRIVATE
    FunctionalParser.cpp
    PairDensitySlater.cpp
    PairDensityVWN.cpp
    PairDensityPBE_C.cpp
    PairDensityPBE_X.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
