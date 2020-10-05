target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NuclearPotentialRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
