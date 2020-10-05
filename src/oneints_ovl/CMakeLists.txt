target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/OverlapIntegralsDriver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForDX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForFG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForGF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForGG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForPX.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OverlapRecFuncForSX.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )