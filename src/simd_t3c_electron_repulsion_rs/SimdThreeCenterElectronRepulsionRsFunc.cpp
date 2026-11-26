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




#include "SimdThreeCenterElectronRepulsionRsFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"

#include "SimdThreeCenterElectronRepulsionRsRecSSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecSIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecPIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecDIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecFIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecGIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecHIL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIPL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIDL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIFL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIGL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHD.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHI.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIHL.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIS.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIP.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIID.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIF.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIG.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIH.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIII.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIK.hpp"
#include "SimdThreeCenterElectronRepulsionRsRecIIL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_electron_repulsion(double               *values,
                              const size_t          npairs,
                              const size_t          natoms,
                              const CBasisFunction &a_function,
                              const CBasisFunction &b_function,
                              const CBasisFunction &c_function,
                              const CSimdMatrix    &coordinates,
                              const CSimdMatrix    &c_coordinates,
                              CSimdMatrix          &buffer,
                              const double          omega,
                              const double          threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    // NOTE: the combination is dispatched on the three angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last. The
    // unattenuated dispatcher does the same and for the same reason.

    switch ((la * 7 + lb) * 9 + lc)
    {
        case     0:
            compute_rs_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     1:
            compute_rs_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     2:
            compute_rs_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     3:
            compute_rs_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     4:
            compute_rs_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     5:
            compute_rs_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     6:
            compute_rs_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     7:
            compute_rs_ssk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     8:
            compute_rs_ssl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case     9:
            compute_rs_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    10:
            compute_rs_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    11:
            compute_rs_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    12:
            compute_rs_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    13:
            compute_rs_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    14:
            compute_rs_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    15:
            compute_rs_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    16:
            compute_rs_spk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    17:
            compute_rs_spl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    18:
            compute_rs_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    19:
            compute_rs_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    20:
            compute_rs_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    21:
            compute_rs_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    22:
            compute_rs_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    23:
            compute_rs_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    24:
            compute_rs_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    25:
            compute_rs_sdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    26:
            compute_rs_sdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    27:
            compute_rs_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    28:
            compute_rs_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    29:
            compute_rs_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    30:
            compute_rs_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    31:
            compute_rs_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    32:
            compute_rs_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    33:
            compute_rs_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    34:
            compute_rs_sfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    35:
            compute_rs_sfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    36:
            compute_rs_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    37:
            compute_rs_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    38:
            compute_rs_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    39:
            compute_rs_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    40:
            compute_rs_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    41:
            compute_rs_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    42:
            compute_rs_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    43:
            compute_rs_sgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    44:
            compute_rs_sgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    45:
            compute_rs_shs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    46:
            compute_rs_shp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    47:
            compute_rs_shd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    48:
            compute_rs_shf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    49:
            compute_rs_shg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    50:
            compute_rs_shh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    51:
            compute_rs_shi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    52:
            compute_rs_shk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    53:
            compute_rs_shl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    54:
            compute_rs_sis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    55:
            compute_rs_sip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    56:
            compute_rs_sid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    57:
            compute_rs_sif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    58:
            compute_rs_sig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    59:
            compute_rs_sih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    60:
            compute_rs_sii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    61:
            compute_rs_sik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    62:
            compute_rs_sil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    63:
            compute_rs_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    64:
            compute_rs_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    65:
            compute_rs_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    66:
            compute_rs_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    67:
            compute_rs_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    68:
            compute_rs_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    69:
            compute_rs_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    70:
            compute_rs_psk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    71:
            compute_rs_psl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    72:
            compute_rs_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    73:
            compute_rs_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    74:
            compute_rs_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    75:
            compute_rs_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    76:
            compute_rs_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    77:
            compute_rs_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    78:
            compute_rs_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    79:
            compute_rs_ppk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    80:
            compute_rs_ppl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    81:
            compute_rs_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    82:
            compute_rs_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    83:
            compute_rs_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    84:
            compute_rs_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    85:
            compute_rs_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    86:
            compute_rs_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    87:
            compute_rs_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    88:
            compute_rs_pdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    89:
            compute_rs_pdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    90:
            compute_rs_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    91:
            compute_rs_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    92:
            compute_rs_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    93:
            compute_rs_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    94:
            compute_rs_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    95:
            compute_rs_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    96:
            compute_rs_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    97:
            compute_rs_pfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    98:
            compute_rs_pfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case    99:
            compute_rs_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   100:
            compute_rs_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   101:
            compute_rs_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   102:
            compute_rs_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   103:
            compute_rs_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   104:
            compute_rs_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   105:
            compute_rs_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   106:
            compute_rs_pgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   107:
            compute_rs_pgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   108:
            compute_rs_phs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   109:
            compute_rs_php_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   110:
            compute_rs_phd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   111:
            compute_rs_phf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   112:
            compute_rs_phg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   113:
            compute_rs_phh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   114:
            compute_rs_phi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   115:
            compute_rs_phk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   116:
            compute_rs_phl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   117:
            compute_rs_pis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   118:
            compute_rs_pip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   119:
            compute_rs_pid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   120:
            compute_rs_pif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   121:
            compute_rs_pig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   122:
            compute_rs_pih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   123:
            compute_rs_pii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   124:
            compute_rs_pik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   125:
            compute_rs_pil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   126:
            compute_rs_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   127:
            compute_rs_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   128:
            compute_rs_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   129:
            compute_rs_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   130:
            compute_rs_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   131:
            compute_rs_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   132:
            compute_rs_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   133:
            compute_rs_dsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   134:
            compute_rs_dsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   135:
            compute_rs_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   136:
            compute_rs_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   137:
            compute_rs_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   138:
            compute_rs_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   139:
            compute_rs_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   140:
            compute_rs_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   141:
            compute_rs_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   142:
            compute_rs_dpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   143:
            compute_rs_dpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   144:
            compute_rs_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   145:
            compute_rs_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   146:
            compute_rs_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   147:
            compute_rs_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   148:
            compute_rs_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   149:
            compute_rs_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   150:
            compute_rs_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   151:
            compute_rs_ddk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   152:
            compute_rs_ddl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   153:
            compute_rs_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   154:
            compute_rs_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   155:
            compute_rs_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   156:
            compute_rs_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   157:
            compute_rs_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   158:
            compute_rs_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   159:
            compute_rs_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   160:
            compute_rs_dfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   161:
            compute_rs_dfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   162:
            compute_rs_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   163:
            compute_rs_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   164:
            compute_rs_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   165:
            compute_rs_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   166:
            compute_rs_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   167:
            compute_rs_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   168:
            compute_rs_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   169:
            compute_rs_dgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   170:
            compute_rs_dgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   171:
            compute_rs_dhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   172:
            compute_rs_dhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   173:
            compute_rs_dhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   174:
            compute_rs_dhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   175:
            compute_rs_dhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   176:
            compute_rs_dhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   177:
            compute_rs_dhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   178:
            compute_rs_dhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   179:
            compute_rs_dhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   180:
            compute_rs_dis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   181:
            compute_rs_dip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   182:
            compute_rs_did_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   183:
            compute_rs_dif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   184:
            compute_rs_dig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   185:
            compute_rs_dih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   186:
            compute_rs_dii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   187:
            compute_rs_dik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   188:
            compute_rs_dil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   189:
            compute_rs_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   190:
            compute_rs_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   191:
            compute_rs_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   192:
            compute_rs_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   193:
            compute_rs_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   194:
            compute_rs_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   195:
            compute_rs_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   196:
            compute_rs_fsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   197:
            compute_rs_fsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   198:
            compute_rs_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   199:
            compute_rs_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   200:
            compute_rs_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   201:
            compute_rs_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   202:
            compute_rs_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   203:
            compute_rs_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   204:
            compute_rs_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   205:
            compute_rs_fpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   206:
            compute_rs_fpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   207:
            compute_rs_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   208:
            compute_rs_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   209:
            compute_rs_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   210:
            compute_rs_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   211:
            compute_rs_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   212:
            compute_rs_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   213:
            compute_rs_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   214:
            compute_rs_fdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   215:
            compute_rs_fdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   216:
            compute_rs_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   217:
            compute_rs_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   218:
            compute_rs_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   219:
            compute_rs_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   220:
            compute_rs_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   221:
            compute_rs_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   222:
            compute_rs_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   223:
            compute_rs_ffk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   224:
            compute_rs_ffl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   225:
            compute_rs_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   226:
            compute_rs_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   227:
            compute_rs_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   228:
            compute_rs_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   229:
            compute_rs_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   230:
            compute_rs_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   231:
            compute_rs_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   232:
            compute_rs_fgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   233:
            compute_rs_fgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   234:
            compute_rs_fhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   235:
            compute_rs_fhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   236:
            compute_rs_fhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   237:
            compute_rs_fhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   238:
            compute_rs_fhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   239:
            compute_rs_fhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   240:
            compute_rs_fhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   241:
            compute_rs_fhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   242:
            compute_rs_fhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   243:
            compute_rs_fis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   244:
            compute_rs_fip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   245:
            compute_rs_fid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   246:
            compute_rs_fif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   247:
            compute_rs_fig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   248:
            compute_rs_fih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   249:
            compute_rs_fii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   250:
            compute_rs_fik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   251:
            compute_rs_fil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   252:
            compute_rs_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   253:
            compute_rs_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   254:
            compute_rs_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   255:
            compute_rs_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   256:
            compute_rs_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   257:
            compute_rs_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   258:
            compute_rs_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   259:
            compute_rs_gsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   260:
            compute_rs_gsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   261:
            compute_rs_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   262:
            compute_rs_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   263:
            compute_rs_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   264:
            compute_rs_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   265:
            compute_rs_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   266:
            compute_rs_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   267:
            compute_rs_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   268:
            compute_rs_gpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   269:
            compute_rs_gpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   270:
            compute_rs_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   271:
            compute_rs_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   272:
            compute_rs_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   273:
            compute_rs_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   274:
            compute_rs_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   275:
            compute_rs_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   276:
            compute_rs_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   277:
            compute_rs_gdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   278:
            compute_rs_gdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   279:
            compute_rs_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   280:
            compute_rs_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   281:
            compute_rs_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   282:
            compute_rs_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   283:
            compute_rs_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   284:
            compute_rs_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   285:
            compute_rs_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   286:
            compute_rs_gfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   287:
            compute_rs_gfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   288:
            compute_rs_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   289:
            compute_rs_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   290:
            compute_rs_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   291:
            compute_rs_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   292:
            compute_rs_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   293:
            compute_rs_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   294:
            compute_rs_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   295:
            compute_rs_ggk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   296:
            compute_rs_ggl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   297:
            compute_rs_ghs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   298:
            compute_rs_ghp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   299:
            compute_rs_ghd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   300:
            compute_rs_ghf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   301:
            compute_rs_ghg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   302:
            compute_rs_ghh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   303:
            compute_rs_ghi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   304:
            compute_rs_ghk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   305:
            compute_rs_ghl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   306:
            compute_rs_gis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   307:
            compute_rs_gip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   308:
            compute_rs_gid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   309:
            compute_rs_gif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   310:
            compute_rs_gig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   311:
            compute_rs_gih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   312:
            compute_rs_gii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   313:
            compute_rs_gik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   314:
            compute_rs_gil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   315:
            compute_rs_hss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   316:
            compute_rs_hsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   317:
            compute_rs_hsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   318:
            compute_rs_hsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   319:
            compute_rs_hsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   320:
            compute_rs_hsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   321:
            compute_rs_hsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   322:
            compute_rs_hsk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   323:
            compute_rs_hsl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   324:
            compute_rs_hps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   325:
            compute_rs_hpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   326:
            compute_rs_hpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   327:
            compute_rs_hpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   328:
            compute_rs_hpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   329:
            compute_rs_hph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   330:
            compute_rs_hpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   331:
            compute_rs_hpk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   332:
            compute_rs_hpl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   333:
            compute_rs_hds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   334:
            compute_rs_hdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   335:
            compute_rs_hdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   336:
            compute_rs_hdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   337:
            compute_rs_hdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   338:
            compute_rs_hdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   339:
            compute_rs_hdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   340:
            compute_rs_hdk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   341:
            compute_rs_hdl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   342:
            compute_rs_hfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   343:
            compute_rs_hfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   344:
            compute_rs_hfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   345:
            compute_rs_hff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   346:
            compute_rs_hfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   347:
            compute_rs_hfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   348:
            compute_rs_hfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   349:
            compute_rs_hfk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   350:
            compute_rs_hfl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   351:
            compute_rs_hgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   352:
            compute_rs_hgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   353:
            compute_rs_hgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   354:
            compute_rs_hgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   355:
            compute_rs_hgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   356:
            compute_rs_hgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   357:
            compute_rs_hgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   358:
            compute_rs_hgk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   359:
            compute_rs_hgl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   360:
            compute_rs_hhs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   361:
            compute_rs_hhp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   362:
            compute_rs_hhd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   363:
            compute_rs_hhf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   364:
            compute_rs_hhg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   365:
            compute_rs_hhh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   366:
            compute_rs_hhi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   367:
            compute_rs_hhk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   368:
            compute_rs_hhl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   369:
            compute_rs_his_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   370:
            compute_rs_hip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   371:
            compute_rs_hid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   372:
            compute_rs_hif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   373:
            compute_rs_hig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   374:
            compute_rs_hih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   375:
            compute_rs_hii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   376:
            compute_rs_hik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   377:
            compute_rs_hil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   378:
            compute_rs_iss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   379:
            compute_rs_isp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   380:
            compute_rs_isd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   381:
            compute_rs_isf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   382:
            compute_rs_isg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   383:
            compute_rs_ish_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   384:
            compute_rs_isi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   385:
            compute_rs_isk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   386:
            compute_rs_isl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   387:
            compute_rs_ips_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   388:
            compute_rs_ipp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   389:
            compute_rs_ipd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   390:
            compute_rs_ipf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   391:
            compute_rs_ipg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   392:
            compute_rs_iph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   393:
            compute_rs_ipi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   394:
            compute_rs_ipk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   395:
            compute_rs_ipl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   396:
            compute_rs_ids_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   397:
            compute_rs_idp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   398:
            compute_rs_idd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   399:
            compute_rs_idf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   400:
            compute_rs_idg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   401:
            compute_rs_idh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   402:
            compute_rs_idi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   403:
            compute_rs_idk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   404:
            compute_rs_idl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   405:
            compute_rs_ifs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   406:
            compute_rs_ifp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   407:
            compute_rs_ifd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   408:
            compute_rs_iff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   409:
            compute_rs_ifg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   410:
            compute_rs_ifh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   411:
            compute_rs_ifi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   412:
            compute_rs_ifk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   413:
            compute_rs_ifl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   414:
            compute_rs_igs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   415:
            compute_rs_igp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   416:
            compute_rs_igd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   417:
            compute_rs_igf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   418:
            compute_rs_igg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   419:
            compute_rs_igh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   420:
            compute_rs_igi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   421:
            compute_rs_igk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   422:
            compute_rs_igl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   423:
            compute_rs_ihs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   424:
            compute_rs_ihp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   425:
            compute_rs_ihd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   426:
            compute_rs_ihf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   427:
            compute_rs_ihg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   428:
            compute_rs_ihh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   429:
            compute_rs_ihi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   430:
            compute_rs_ihk_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   431:
            compute_rs_ihl_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   432:
            compute_rs_iis_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   433:
            compute_rs_iip_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   434:
            compute_rs_iid_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   435:
            compute_rs_iif_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   436:
            compute_rs_iig_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   437:
            compute_rs_iih_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   438:
            compute_rs_iii_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   439:
            compute_rs_iik_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        case   440:
            compute_rs_iil_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates, c_coordinates,
                buffer, omega, threshold);
            return;

        default:
            break;
    }

    // NOTE: the combination has no kernel. The dispatcher stops rather than leaving
    // the values of the tensor unwritten, which is what a caller would otherwise read
    // as integrals. The range separated kernels reach angular momentum six on the two
    // orbital sides and eight on the auxiliary one, as the unattenuated ones do.

    errors::assertMsgCritical(
        false,
        std::string("SimdThreeCenterElectronRepulsionRsFunc.compute_rs_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simdt3ceri
