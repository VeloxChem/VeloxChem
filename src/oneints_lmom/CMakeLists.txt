target_sources(vlxobjs
  PRIVATE
    LinearMomentumIntegralsDriver.cpp
    LinearMomentumRecFuncForDX.cpp
    LinearMomentumRecFuncForFF.cpp
    LinearMomentumRecFuncForFG.cpp
    LinearMomentumRecFuncForGF.cpp
    LinearMomentumRecFuncForGG.cpp
    LinearMomentumRecFuncForPX.cpp
    LinearMomentumRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
