target_sources(vlxobjs
  PRIVATE
    AngularMomentumIntegrals.cpp
    BoysFuncTable.cpp
    ElectricFieldIntegrals.cpp
    ElectricFieldValues.cpp
    ElectricFieldGradients.cpp
    GtoInfo.cpp
    LinearMomentumIntegrals.cpp
    NuclearPotentialValues.cpp
    QuadrupoleIntegrals.cpp
    OldOneElecIntsDrivers.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
