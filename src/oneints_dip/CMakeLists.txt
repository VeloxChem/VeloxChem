target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricDipoleRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
