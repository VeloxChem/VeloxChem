target_sources(vlxobjs
  PRIVATE
    ElectronRepulsionGeom1000ContrRecDDXX.cpp
    ElectronRepulsionGeom1000ContrRecDFXX.cpp
    ElectronRepulsionGeom1000ContrRecDGXX.cpp
    ElectronRepulsionGeom1000ContrRecDHXX.cpp
    ElectronRepulsionGeom1000ContrRecDIXX.cpp
    ElectronRepulsionGeom1000ContrRecDPXX.cpp
    ElectronRepulsionGeom1000ContrRecDSXX.cpp
    ElectronRepulsionGeom1000ContrRecFDXX.cpp
    ElectronRepulsionGeom1000ContrRecFFXX.cpp
    ElectronRepulsionGeom1000ContrRecFGXX.cpp
    ElectronRepulsionGeom1000ContrRecFHXX.cpp
    ElectronRepulsionGeom1000ContrRecFIXX.cpp
    ElectronRepulsionGeom1000ContrRecFPXX.cpp
    ElectronRepulsionGeom1000ContrRecFSXX.cpp
    ElectronRepulsionGeom1000ContrRecGDXX.cpp
    ElectronRepulsionGeom1000ContrRecGFXX.cpp
    ElectronRepulsionGeom1000ContrRecGGXX.cpp
    ElectronRepulsionGeom1000ContrRecGHXX.cpp
    ElectronRepulsionGeom1000ContrRecGIXX.cpp
    ElectronRepulsionGeom1000ContrRecGPXX.cpp
    ElectronRepulsionGeom1000ContrRecGSXX.cpp
    ElectronRepulsionGeom1000ContrRecHDXX.cpp
    ElectronRepulsionGeom1000ContrRecHFXX.cpp
    ElectronRepulsionGeom1000ContrRecHGXX.cpp
    ElectronRepulsionGeom1000ContrRecHHXX.cpp
    ElectronRepulsionGeom1000ContrRecHIXX.cpp
    ElectronRepulsionGeom1000ContrRecHPXX.cpp
    ElectronRepulsionGeom1000ContrRecHSXX.cpp
    ElectronRepulsionGeom1000ContrRecIDXX.cpp
    ElectronRepulsionGeom1000ContrRecIFXX.cpp
    ElectronRepulsionGeom1000ContrRecIGXX.cpp
    ElectronRepulsionGeom1000ContrRecIHXX.cpp
    ElectronRepulsionGeom1000ContrRecIIXX.cpp
    ElectronRepulsionGeom1000ContrRecIPXX.cpp
    ElectronRepulsionGeom1000ContrRecISXX.cpp
    ElectronRepulsionGeom1000ContrRecPDXX.cpp
    ElectronRepulsionGeom1000ContrRecPFXX.cpp
    ElectronRepulsionGeom1000ContrRecPGXX.cpp
    ElectronRepulsionGeom1000ContrRecPHXX.cpp
    ElectronRepulsionGeom1000ContrRecPIXX.cpp
    ElectronRepulsionGeom1000ContrRecPPXX.cpp
    ElectronRepulsionGeom1000ContrRecPSXX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
