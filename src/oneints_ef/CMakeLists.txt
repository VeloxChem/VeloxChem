target_sources(vlxobjs
  PRIVATE
    ElectricFieldIntegralsDriver.cpp
    ElectricFieldRecFuncForDX.cpp
    ElectricFieldRecFuncForFF.cpp
    ElectricFieldRecFuncForFG.cpp
    ElectricFieldRecFuncForGF.cpp
    ElectricFieldRecFuncForGG.cpp
    ElectricFieldRecFuncForPX.cpp
    ElectricFieldRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
