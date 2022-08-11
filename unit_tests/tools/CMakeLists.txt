target_sources(utests
  PRIVATE
    CheckFunctions.cpp
    DummyFunctionals.cpp
    DummyFunctions.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
