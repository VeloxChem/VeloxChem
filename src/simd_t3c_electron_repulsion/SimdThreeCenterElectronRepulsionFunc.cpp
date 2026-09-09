//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#include "SimdThreeCenterElectronRepulsionFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDID.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDII.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFID.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFII.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGID.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGII.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHID.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHII.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIID.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIII.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecIPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPID.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPII.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_electron_repulsion(double               *values,
                           const size_t          npairs,
                           const size_t          natoms,
                           const CBasisFunction &a_function,
                           const CBasisFunction &b_function,
                           const CBasisFunction &c_function,
                           const CSimdMatrix    &coordinates,
                           const CSimdMatrix    &c_coordinates,
                           CSimdMatrix          &buffer,
                           const double          threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    // NOTE: the combination is dispatched on the three angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last.

    switch ((la * 7 + lb) * 9 + lc)
    {
        case     0:
            compute_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     1:
            compute_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     2:
            compute_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     3:
            compute_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     4:
            compute_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     5:
            compute_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     6:
            compute_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     7:
            compute_ssk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     8:
            compute_ssl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case     9:
            compute_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    10:
            compute_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    11:
            compute_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    12:
            compute_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    13:
            compute_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    14:
            compute_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    15:
            compute_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    16:
            compute_spk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    17:
            compute_spl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    18:
            compute_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    19:
            compute_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    20:
            compute_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    21:
            compute_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    22:
            compute_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    23:
            compute_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    24:
            compute_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    25:
            compute_sdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    26:
            compute_sdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    27:
            compute_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    28:
            compute_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    29:
            compute_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    30:
            compute_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    31:
            compute_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    32:
            compute_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    33:
            compute_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    34:
            compute_sfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    35:
            compute_sfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    36:
            compute_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    37:
            compute_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    38:
            compute_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    39:
            compute_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    40:
            compute_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    41:
            compute_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    42:
            compute_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    43:
            compute_sgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    44:
            compute_sgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    45:
            compute_shs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    46:
            compute_shp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    47:
            compute_shd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    48:
            compute_shf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    49:
            compute_shg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    50:
            compute_shh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    51:
            compute_shi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    52:
            compute_shk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    53:
            compute_shl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    54:
            compute_sis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    55:
            compute_sip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    56:
            compute_sid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    57:
            compute_sif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    58:
            compute_sig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    59:
            compute_sih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    60:
            compute_sii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    61:
            compute_sik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    62:
            compute_sil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    63:
            compute_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    64:
            compute_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    65:
            compute_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    66:
            compute_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    67:
            compute_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    68:
            compute_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    69:
            compute_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    70:
            compute_psk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    71:
            compute_psl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    72:
            compute_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    73:
            compute_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    74:
            compute_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    75:
            compute_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    76:
            compute_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    77:
            compute_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    78:
            compute_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    79:
            compute_ppk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    80:
            compute_ppl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    81:
            compute_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    82:
            compute_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    83:
            compute_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    84:
            compute_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    85:
            compute_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    86:
            compute_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    87:
            compute_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    88:
            compute_pdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    89:
            compute_pdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    90:
            compute_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    91:
            compute_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    92:
            compute_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    93:
            compute_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    94:
            compute_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    95:
            compute_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    96:
            compute_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    97:
            compute_pfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    98:
            compute_pfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case    99:
            compute_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   100:
            compute_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   101:
            compute_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   102:
            compute_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   103:
            compute_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   104:
            compute_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   105:
            compute_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   106:
            compute_pgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   107:
            compute_pgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   108:
            compute_phs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   109:
            compute_php_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   110:
            compute_phd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   111:
            compute_phf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   112:
            compute_phg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   113:
            compute_phh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   114:
            compute_phi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   115:
            compute_phk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   116:
            compute_phl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   117:
            compute_pis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   118:
            compute_pip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   119:
            compute_pid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   120:
            compute_pif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   121:
            compute_pig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   122:
            compute_pih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   123:
            compute_pii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   124:
            compute_pik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   125:
            compute_pil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   126:
            compute_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   127:
            compute_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   128:
            compute_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   129:
            compute_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   130:
            compute_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   131:
            compute_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   132:
            compute_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   133:
            compute_dsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   134:
            compute_dsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   135:
            compute_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   136:
            compute_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   137:
            compute_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   138:
            compute_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   139:
            compute_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   140:
            compute_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   141:
            compute_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   142:
            compute_dpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   143:
            compute_dpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   144:
            compute_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   145:
            compute_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   146:
            compute_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   147:
            compute_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   148:
            compute_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   149:
            compute_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   150:
            compute_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   151:
            compute_ddk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   152:
            compute_ddl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   153:
            compute_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   154:
            compute_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   155:
            compute_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   156:
            compute_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   157:
            compute_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   158:
            compute_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   159:
            compute_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   160:
            compute_dfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   161:
            compute_dfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   162:
            compute_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   163:
            compute_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   164:
            compute_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   165:
            compute_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   166:
            compute_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   167:
            compute_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   168:
            compute_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   169:
            compute_dgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   170:
            compute_dgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   171:
            compute_dhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   172:
            compute_dhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   173:
            compute_dhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   174:
            compute_dhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   175:
            compute_dhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   176:
            compute_dhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   177:
            compute_dhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   178:
            compute_dhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   179:
            compute_dhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   180:
            compute_dis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   181:
            compute_dip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   182:
            compute_did_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   183:
            compute_dif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   184:
            compute_dig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   185:
            compute_dih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   186:
            compute_dii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   187:
            compute_dik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   188:
            compute_dil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   189:
            compute_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   190:
            compute_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   191:
            compute_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   192:
            compute_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   193:
            compute_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   194:
            compute_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   195:
            compute_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   196:
            compute_fsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   197:
            compute_fsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   198:
            compute_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   199:
            compute_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   200:
            compute_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   201:
            compute_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   202:
            compute_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   203:
            compute_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   204:
            compute_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   205:
            compute_fpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   206:
            compute_fpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   207:
            compute_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   208:
            compute_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   209:
            compute_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   210:
            compute_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   211:
            compute_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   212:
            compute_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   213:
            compute_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   214:
            compute_fdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   215:
            compute_fdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   216:
            compute_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   217:
            compute_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   218:
            compute_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   219:
            compute_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   220:
            compute_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   221:
            compute_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   222:
            compute_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   223:
            compute_ffk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   224:
            compute_ffl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   225:
            compute_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   226:
            compute_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   227:
            compute_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   228:
            compute_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   229:
            compute_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   230:
            compute_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   231:
            compute_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   232:
            compute_fgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   233:
            compute_fgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   234:
            compute_fhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   235:
            compute_fhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   236:
            compute_fhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   237:
            compute_fhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   238:
            compute_fhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   239:
            compute_fhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   240:
            compute_fhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   241:
            compute_fhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   242:
            compute_fhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   243:
            compute_fis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   244:
            compute_fip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   245:
            compute_fid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   246:
            compute_fif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   247:
            compute_fig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   248:
            compute_fih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   249:
            compute_fii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   250:
            compute_fik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   251:
            compute_fil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   252:
            compute_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   253:
            compute_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   254:
            compute_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   255:
            compute_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   256:
            compute_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   257:
            compute_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   258:
            compute_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   259:
            compute_gsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   260:
            compute_gsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   261:
            compute_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   262:
            compute_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   263:
            compute_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   264:
            compute_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   265:
            compute_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   266:
            compute_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   267:
            compute_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   268:
            compute_gpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   269:
            compute_gpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   270:
            compute_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   271:
            compute_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   272:
            compute_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   273:
            compute_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   274:
            compute_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   275:
            compute_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   276:
            compute_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   277:
            compute_gdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   278:
            compute_gdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   279:
            compute_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   280:
            compute_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   281:
            compute_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   282:
            compute_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   283:
            compute_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   284:
            compute_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   285:
            compute_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   286:
            compute_gfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   287:
            compute_gfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   288:
            compute_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   289:
            compute_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   290:
            compute_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   291:
            compute_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   292:
            compute_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   293:
            compute_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   294:
            compute_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   295:
            compute_ggk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   296:
            compute_ggl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   297:
            compute_ghs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   298:
            compute_ghp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   299:
            compute_ghd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   300:
            compute_ghf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   301:
            compute_ghg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   302:
            compute_ghh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   303:
            compute_ghi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   304:
            compute_ghk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   305:
            compute_ghl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   306:
            compute_gis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   307:
            compute_gip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   308:
            compute_gid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   309:
            compute_gif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   310:
            compute_gig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   311:
            compute_gih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   312:
            compute_gii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   313:
            compute_gik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   314:
            compute_gil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   315:
            compute_hss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   316:
            compute_hsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   317:
            compute_hsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   318:
            compute_hsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   319:
            compute_hsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   320:
            compute_hsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   321:
            compute_hsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   322:
            compute_hsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   323:
            compute_hsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   324:
            compute_hps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   325:
            compute_hpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   326:
            compute_hpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   327:
            compute_hpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   328:
            compute_hpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   329:
            compute_hph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   330:
            compute_hpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   331:
            compute_hpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   332:
            compute_hpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   333:
            compute_hds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   334:
            compute_hdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   335:
            compute_hdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   336:
            compute_hdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   337:
            compute_hdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   338:
            compute_hdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   339:
            compute_hdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   340:
            compute_hdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   341:
            compute_hdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   342:
            compute_hfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   343:
            compute_hfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   344:
            compute_hfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   345:
            compute_hff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   346:
            compute_hfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   347:
            compute_hfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   348:
            compute_hfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   349:
            compute_hfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   350:
            compute_hfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   351:
            compute_hgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   352:
            compute_hgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   353:
            compute_hgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   354:
            compute_hgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   355:
            compute_hgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   356:
            compute_hgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   357:
            compute_hgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   358:
            compute_hgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   359:
            compute_hgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   360:
            compute_hhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   361:
            compute_hhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   362:
            compute_hhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   363:
            compute_hhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   364:
            compute_hhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   365:
            compute_hhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   366:
            compute_hhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   367:
            compute_hhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   368:
            compute_hhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   369:
            compute_his_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   370:
            compute_hip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   371:
            compute_hid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   372:
            compute_hif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   373:
            compute_hig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   374:
            compute_hih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   375:
            compute_hii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   376:
            compute_hik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   377:
            compute_hil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   378:
            compute_iss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   379:
            compute_isp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   380:
            compute_isd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   381:
            compute_isf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   382:
            compute_isg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   383:
            compute_ish_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   384:
            compute_isi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   385:
            compute_isk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   386:
            compute_isl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   387:
            compute_ips_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   388:
            compute_ipp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   389:
            compute_ipd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   390:
            compute_ipf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   391:
            compute_ipg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   392:
            compute_iph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   393:
            compute_ipi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   394:
            compute_ipk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   395:
            compute_ipl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   396:
            compute_ids_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   397:
            compute_idp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   398:
            compute_idd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   399:
            compute_idf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   400:
            compute_idg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   401:
            compute_idh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   402:
            compute_idi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   403:
            compute_idk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   404:
            compute_idl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   405:
            compute_ifs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   406:
            compute_ifp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   407:
            compute_ifd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   408:
            compute_iff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   409:
            compute_ifg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   410:
            compute_ifh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   411:
            compute_ifi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   412:
            compute_ifk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   413:
            compute_ifl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   414:
            compute_igs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   415:
            compute_igp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   416:
            compute_igd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   417:
            compute_igf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   418:
            compute_igg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   419:
            compute_igh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   420:
            compute_igi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   421:
            compute_igk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   422:
            compute_igl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   423:
            compute_ihs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   424:
            compute_ihp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   425:
            compute_ihd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   426:
            compute_ihf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   427:
            compute_ihg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   428:
            compute_ihh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   429:
            compute_ihi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   430:
            compute_ihk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   431:
            compute_ihl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   432:
            compute_iis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   433:
            compute_iip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   434:
            compute_iid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   435:
            compute_iif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   436:
            compute_iig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   437:
            compute_iih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   438:
            compute_iii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   439:
            compute_iik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        case   440:
            compute_iil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates, buffer, threshold);
            return;

        default:
            break;
    }

    // NOTE: the combination has no kernel. The dispatcher stops rather than leaving
    // the values of the tensor unwritten, which is what a caller would otherwise read
    // as integrals.

    errors::assertMsgCritical(
        false, std::string("SimdThreeCenterElectronRepulsionFunc.compute_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simdt3ceri
