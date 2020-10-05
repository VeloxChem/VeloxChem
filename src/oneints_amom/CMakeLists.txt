target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AngularMomentumRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
