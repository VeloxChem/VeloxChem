target_sources(vlxobjs
  PRIVATE
    ElectronRepulsionGeom0010ContrRecXXDD.cpp
    ElectronRepulsionGeom0010ContrRecXXDP.cpp
    ElectronRepulsionGeom0010ContrRecXXDS.cpp
    ElectronRepulsionGeom0010ContrRecXXPD.cpp
    ElectronRepulsionGeom0010ContrRecXXPF.cpp
    ElectronRepulsionGeom0010ContrRecXXPP.cpp
    ElectronRepulsionGeom0010ContrRecXXPS.cpp
    ElectronRepulsionGeom0010ContrRecXXSD.cpp
    ElectronRepulsionGeom0010ContrRecXXSF.cpp
    ElectronRepulsionGeom0010ContrRecXXSG.cpp
    ElectronRepulsionGeom0010ContrRecXXSP.cpp
    ElectronRepulsionGeom0010ContrRecXXSS.cpp
    ElectronRepulsionGeom0100ContrRecDDXX.cpp
    ElectronRepulsionGeom0100ContrRecDPXX.cpp
    ElectronRepulsionGeom0100ContrRecPDXX.cpp
    ElectronRepulsionGeom0100ContrRecPFXX.cpp
    ElectronRepulsionGeom0100ContrRecPPXX.cpp
    ElectronRepulsionGeom0100ContrRecPSXX.cpp
    ElectronRepulsionGeom0100ContrRecSDXX.cpp
    ElectronRepulsionGeom0100ContrRecSFXX.cpp
    ElectronRepulsionGeom0100ContrRecSGXX.cpp
    ElectronRepulsionGeom0100ContrRecSHXX.cpp
    ElectronRepulsionGeom0100ContrRecSPXX.cpp
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
    ElectronRepulsionGeom1010ContrRecDDXX.cpp
    ElectronRepulsionGeom1010ContrRecDPXX.cpp
    ElectronRepulsionGeom1010ContrRecDSXX.cpp
    ElectronRepulsionGeom1010ContrRecPDXX.cpp
    ElectronRepulsionGeom1010ContrRecPFXX.cpp
    ElectronRepulsionGeom1010ContrRecPPXX.cpp
    ElectronRepulsionGeom1010ContrRecPSXX.cpp
    ElectronRepulsionGeom1010ContrRecSDXX.cpp
    ElectronRepulsionGeom1010ContrRecSFXX.cpp
    ElectronRepulsionGeom1010ContrRecSGXX.cpp
    ElectronRepulsionGeom1010ContrRecSPXX.cpp
    ElectronRepulsionGeom1010ContrRecSSXX.cpp
    ElectronRepulsionGeom1100ContrRecDDXX.cpp
    ElectronRepulsionGeom1100ContrRecDPXX.cpp
    ElectronRepulsionGeom1100ContrRecDSXX.cpp
    ElectronRepulsionGeom1100ContrRecPDXX.cpp
    ElectronRepulsionGeom1100ContrRecPFXX.cpp
    ElectronRepulsionGeom1100ContrRecPPXX.cpp
    ElectronRepulsionGeom1100ContrRecPSXX.cpp
    ElectronRepulsionGeom1100ContrRecSDXX.cpp
    ElectronRepulsionGeom1100ContrRecSFXX.cpp
    ElectronRepulsionGeom1100ContrRecSGXX.cpp
    ElectronRepulsionGeom1100ContrRecSPXX.cpp
    ElectronRepulsionGeom2000ContrRecDDXX.cpp
    ElectronRepulsionGeom2000ContrRecDPXX.cpp
    ElectronRepulsionGeom2000ContrRecDSXX.cpp
    ElectronRepulsionGeom2000ContrRecPDXX.cpp
    ElectronRepulsionGeom2000ContrRecPFXX.cpp
    ElectronRepulsionGeom2000ContrRecPPXX.cpp
    ElectronRepulsionGeom2000ContrRecPSXX.cpp
    ElectronRepulsionGeom2000ContrRecSXXX.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
