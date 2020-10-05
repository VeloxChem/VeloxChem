target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/KineticEnergyRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
