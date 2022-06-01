target_sources(vlxobjs
  PRIVATE
    KineticEnergyIntegralsDriver.cpp
    KineticEnergyRecFuncForDX.cpp
    KineticEnergyRecFuncForFF.cpp
    KineticEnergyRecFuncForFG.cpp
    KineticEnergyRecFuncForGF.cpp
    KineticEnergyRecFuncForGG.cpp
    KineticEnergyRecFuncForPX.cpp
    KineticEnergyRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
