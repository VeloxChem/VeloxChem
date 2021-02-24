target_sources(vlxobjs
  PRIVATE
    ElectricDipoleIntegralsDriver.cpp
    ElectricDipoleRecFuncForDX.cpp
    ElectricDipoleRecFuncForFF.cpp
    ElectricDipoleRecFuncForFG.cpp
    ElectricDipoleRecFuncForGF.cpp
    ElectricDipoleRecFuncForGG.cpp
    ElectricDipoleRecFuncForPX.cpp
    ElectricDipoleRecFuncForSX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
