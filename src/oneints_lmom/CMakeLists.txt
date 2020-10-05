target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LinearMomentumRecFuncForSX.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )