target_sources(vlxobjs
  PRIVATE
    Becke88Functional.cpp
    FunctionalParser.cpp
    NewFunctionalParser.cpp
    LYPFunctional.cpp
    PkzbFunctional.cpp
    SlaterFunctional.cpp
    VWN3Functional.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
