target_sources(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/CheckFunctions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/DummyFunctionals.cpp
    ${CMAKE_CURRENT_LIST_DIR}/DummyFunctions.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
