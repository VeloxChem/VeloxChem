target_sources(vlxobjs
  PRIVATE
    ThreeCenterElectronRepulsionGeom100ContrRecDXX.cpp
    ThreeCenterElectronRepulsionGeom100ContrRecFXX.cpp
    ThreeCenterElectronRepulsionGeom100ContrRecGXX.cpp
    ThreeCenterElectronRepulsionGeom100ContrRecHXX.cpp
    ThreeCenterElectronRepulsionGeom100ContrRecIXX.cpp
    ThreeCenterElectronRepulsionGeom100ContrRecPXX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
