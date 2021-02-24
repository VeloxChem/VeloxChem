target_sources(vlxobjs
  PRIVATE
    NuclearPotentialIntegralsDriver.cpp
    NuclearPotentialRecFuncForDX.cpp
    NuclearPotentialRecFuncForFF.cpp
    NuclearPotentialRecFuncForFG.cpp
    NuclearPotentialRecFuncForGF.cpp
    NuclearPotentialRecFuncForGG.cpp
    NuclearPotentialRecFuncForPX.cpp
    NuclearPotentialRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
