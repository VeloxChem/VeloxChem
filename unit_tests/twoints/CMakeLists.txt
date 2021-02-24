target_sources(utests
  PRIVATE
    AODensityMatrixSetter.cpp
    AOFockMatrixTest.cpp
    CauchySchwarzScreenerTest.cpp
    ElectronRepulsionIntegralsDriverTest.cpp
    FockContainerTest.cpp
    FockSubMatrixTest.cpp
    MOIntsBatchTest.cpp
    ScreeningContainerTest.cpp
    TwoIntsDistributorTest.cpp
    TwoIntsFuncTest.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
