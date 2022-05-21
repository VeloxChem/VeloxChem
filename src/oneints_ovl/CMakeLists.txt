target_sources(vlxobjs
  PRIVATE
    OverlapIntegralsDriver.cpp
    OverlapRecFuncForDX.cpp
    OverlapRecFuncForFF.cpp
    OverlapRecFuncForFG.cpp
    OverlapRecFuncForGF.cpp
    OverlapRecFuncForGG.cpp
    OverlapRecFuncForPX.cpp
    OverlapRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
