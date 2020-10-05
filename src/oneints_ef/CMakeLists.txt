target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectricFieldRecFuncForSX.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )