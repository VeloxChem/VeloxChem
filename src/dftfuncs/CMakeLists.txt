target_sources(vlxobjs
  PRIVATE
    NewFunctionalParser.cpp
    PairDensitySlater.cpp
    PairDensityVWN.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
