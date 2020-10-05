target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/BoysFunction.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )