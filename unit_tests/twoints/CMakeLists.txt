target_sources(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/AODensityMatrixSetter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AOFockMatrixTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CauchySchwarzScreenerTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ElectronRepulsionIntegralsDriverTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FockContainerTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FockSubMatrixTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MOIntsBatchTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScreeningContainerTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TwoIntsDistributorTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TwoIntsFuncTest.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
