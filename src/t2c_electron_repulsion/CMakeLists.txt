target_sources(vlxobjs
  PRIVATE
    TwoCenterElectronRepulsionDriver.cpp
    TwoCenterElectronRepulsionPrimRecDD.cpp
    TwoCenterElectronRepulsionPrimRecDF.cpp
    TwoCenterElectronRepulsionPrimRecDG.cpp
    TwoCenterElectronRepulsionPrimRecDH.cpp
    TwoCenterElectronRepulsionPrimRecDI.cpp
    TwoCenterElectronRepulsionPrimRecDP.cpp
    TwoCenterElectronRepulsionPrimRecDS.cpp
    TwoCenterElectronRepulsionPrimRecFD.cpp
    TwoCenterElectronRepulsionPrimRecFF.cpp
    TwoCenterElectronRepulsionPrimRecFG.cpp
    TwoCenterElectronRepulsionPrimRecFH.cpp
    TwoCenterElectronRepulsionPrimRecFI.cpp
    TwoCenterElectronRepulsionPrimRecFP.cpp
    TwoCenterElectronRepulsionPrimRecFS.cpp
    TwoCenterElectronRepulsionPrimRecGD.cpp
    TwoCenterElectronRepulsionPrimRecGF.cpp
    TwoCenterElectronRepulsionPrimRecGG.cpp
    TwoCenterElectronRepulsionPrimRecGH.cpp
    TwoCenterElectronRepulsionPrimRecGI.cpp
    TwoCenterElectronRepulsionPrimRecGP.cpp
    TwoCenterElectronRepulsionPrimRecGS.cpp
    TwoCenterElectronRepulsionPrimRecHD.cpp
    TwoCenterElectronRepulsionPrimRecHF.cpp
    TwoCenterElectronRepulsionPrimRecHG.cpp
    TwoCenterElectronRepulsionPrimRecHH.cpp
    TwoCenterElectronRepulsionPrimRecHI.cpp
    TwoCenterElectronRepulsionPrimRecHP.cpp
    TwoCenterElectronRepulsionPrimRecHS.cpp
    TwoCenterElectronRepulsionPrimRecID.cpp
    TwoCenterElectronRepulsionPrimRecIF.cpp
    TwoCenterElectronRepulsionPrimRecIG.cpp
    TwoCenterElectronRepulsionPrimRecIH.cpp
    TwoCenterElectronRepulsionPrimRecII.cpp
    TwoCenterElectronRepulsionPrimRecIP.cpp
    TwoCenterElectronRepulsionPrimRecIS.cpp
    TwoCenterElectronRepulsionPrimRecKD.cpp
    TwoCenterElectronRepulsionPrimRecKF.cpp
    TwoCenterElectronRepulsionPrimRecKG.cpp
    TwoCenterElectronRepulsionPrimRecKH.cpp
    TwoCenterElectronRepulsionPrimRecKI.cpp
    TwoCenterElectronRepulsionPrimRecKP.cpp
    TwoCenterElectronRepulsionPrimRecKS.cpp
    TwoCenterElectronRepulsionPrimRecPD.cpp
    TwoCenterElectronRepulsionPrimRecPF.cpp
    TwoCenterElectronRepulsionPrimRecPG.cpp
    TwoCenterElectronRepulsionPrimRecPH.cpp
    TwoCenterElectronRepulsionPrimRecPI.cpp
    TwoCenterElectronRepulsionPrimRecPP.cpp
    TwoCenterElectronRepulsionPrimRecPS.cpp
    TwoCenterElectronRepulsionPrimRecSD.cpp
    TwoCenterElectronRepulsionPrimRecSF.cpp
    TwoCenterElectronRepulsionPrimRecSG.cpp
    TwoCenterElectronRepulsionPrimRecSH.cpp
    TwoCenterElectronRepulsionPrimRecSI.cpp
    TwoCenterElectronRepulsionPrimRecSP.cpp
    TwoCenterElectronRepulsionPrimRecSS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
