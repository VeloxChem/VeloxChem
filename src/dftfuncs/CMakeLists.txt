target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Becke88Functional.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FunctionalParser.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LYPFunctional.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SlaterFunctional.cpp
    ${CMAKE_CURRENT_LIST_DIR}/VWN3Functional.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )