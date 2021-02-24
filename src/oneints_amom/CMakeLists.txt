target_sources(vlxobjs
  PRIVATE
    AngularMomentumIntegralsDriver.cpp
    AngularMomentumRecFuncForDX.cpp
    AngularMomentumRecFuncForFF.cpp
    AngularMomentumRecFuncForFG.cpp
    AngularMomentumRecFuncForGF.cpp
    AngularMomentumRecFuncForGG.cpp
    AngularMomentumRecFuncForPX.cpp
    AngularMomentumRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
