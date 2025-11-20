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

#ifndef ElectronRepulsionGeom1010Func_hpp
#define ElectronRepulsionGeom1010Func_hpp

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom1010RecSSSS.hpp"
#include "ElectronRepulsionGeom1010RecSSSP.hpp"
#include "ElectronRepulsionGeom1010RecSSPS.hpp"
#include "ElectronRepulsionGeom1010RecSSSD.hpp"
#include "ElectronRepulsionGeom1010RecSSDS.hpp"
#include "ElectronRepulsionGeom1010RecSSSF.hpp"
#include "ElectronRepulsionGeom1010RecSSFS.hpp"
#include "ElectronRepulsionGeom1010RecSSPP.hpp"
#include "ElectronRepulsionGeom1010RecSSPD.hpp"
#include "ElectronRepulsionGeom1010RecSSDP.hpp"
#include "ElectronRepulsionGeom1010RecSSPF.hpp"
#include "ElectronRepulsionGeom1010RecSSFP.hpp"
#include "ElectronRepulsionGeom1010RecSSDD.hpp"
#include "ElectronRepulsionGeom1010RecSSDF.hpp"
#include "ElectronRepulsionGeom1010RecSSFD.hpp"
#include "ElectronRepulsionGeom1010RecSSFF.hpp"

#include "ElectronRepulsionGeom1010RecSPSS.hpp"
#include "ElectronRepulsionGeom1010RecSPSP.hpp"
#include "ElectronRepulsionGeom1010RecSPPS.hpp"
#include "ElectronRepulsionGeom1010RecSPSD.hpp"
#include "ElectronRepulsionGeom1010RecSPDS.hpp"
#include "ElectronRepulsionGeom1010RecSPSF.hpp"
#include "ElectronRepulsionGeom1010RecSPFS.hpp"
#include "ElectronRepulsionGeom1010RecSPPP.hpp"
#include "ElectronRepulsionGeom1010RecSPPD.hpp"
#include "ElectronRepulsionGeom1010RecSPDP.hpp"
#include "ElectronRepulsionGeom1010RecSPPF.hpp"
#include "ElectronRepulsionGeom1010RecSPFP.hpp"
#include "ElectronRepulsionGeom1010RecSPDD.hpp"
#include "ElectronRepulsionGeom1010RecSPDF.hpp"
#include "ElectronRepulsionGeom1010RecSPFD.hpp"
#include "ElectronRepulsionGeom1010RecSPFF.hpp"

#include "ElectronRepulsionGeom1010RecPSSS.hpp"
#include "ElectronRepulsionGeom1010RecPSSP.hpp"
#include "ElectronRepulsionGeom1010RecPSPS.hpp"
#include "ElectronRepulsionGeom1010RecPSSD.hpp"
#include "ElectronRepulsionGeom1010RecPSDS.hpp"
#include "ElectronRepulsionGeom1010RecPSSF.hpp"
#include "ElectronRepulsionGeom1010RecPSFS.hpp"
#include "ElectronRepulsionGeom1010RecPSPP.hpp"
#include "ElectronRepulsionGeom1010RecPSPD.hpp"
#include "ElectronRepulsionGeom1010RecPSDP.hpp"
#include "ElectronRepulsionGeom1010RecPSPF.hpp"
#include "ElectronRepulsionGeom1010RecPSFP.hpp"
#include "ElectronRepulsionGeom1010RecPSDD.hpp"
#include "ElectronRepulsionGeom1010RecPSDF.hpp"
#include "ElectronRepulsionGeom1010RecPSFD.hpp"
#include "ElectronRepulsionGeom1010RecPSFF.hpp"

#include "ElectronRepulsionGeom1010RecSDSS.hpp"
#include "ElectronRepulsionGeom1010RecSDSP.hpp"
#include "ElectronRepulsionGeom1010RecSDPS.hpp"
#include "ElectronRepulsionGeom1010RecSDSD.hpp"
#include "ElectronRepulsionGeom1010RecSDDS.hpp"
#include "ElectronRepulsionGeom1010RecSDSF.hpp"
#include "ElectronRepulsionGeom1010RecSDFS.hpp"
#include "ElectronRepulsionGeom1010RecSDPP.hpp"
#include "ElectronRepulsionGeom1010RecSDPD.hpp"
#include "ElectronRepulsionGeom1010RecSDDP.hpp"
#include "ElectronRepulsionGeom1010RecSDPF.hpp"
#include "ElectronRepulsionGeom1010RecSDFP.hpp"
#include "ElectronRepulsionGeom1010RecSDDD.hpp"
#include "ElectronRepulsionGeom1010RecSDDF.hpp"
#include "ElectronRepulsionGeom1010RecSDFD.hpp"
#include "ElectronRepulsionGeom1010RecSDFF.hpp"

#include "ElectronRepulsionGeom1010RecDSSS.hpp"
#include "ElectronRepulsionGeom1010RecDSSP.hpp"
#include "ElectronRepulsionGeom1010RecDSPS.hpp"
#include "ElectronRepulsionGeom1010RecDSSD.hpp"
#include "ElectronRepulsionGeom1010RecDSDS.hpp"
#include "ElectronRepulsionGeom1010RecDSSF.hpp"
#include "ElectronRepulsionGeom1010RecDSFS.hpp"
#include "ElectronRepulsionGeom1010RecDSPP.hpp"
#include "ElectronRepulsionGeom1010RecDSPD.hpp"
#include "ElectronRepulsionGeom1010RecDSDP.hpp"
#include "ElectronRepulsionGeom1010RecDSPF.hpp"
#include "ElectronRepulsionGeom1010RecDSFP.hpp"
#include "ElectronRepulsionGeom1010RecDSDD.hpp"
#include "ElectronRepulsionGeom1010RecDSDF.hpp"
#include "ElectronRepulsionGeom1010RecDSFD.hpp"
#include "ElectronRepulsionGeom1010RecDSFF.hpp"

#include "ElectronRepulsionGeom1010RecSFSS.hpp"
#include "ElectronRepulsionGeom1010RecSFSP.hpp"
#include "ElectronRepulsionGeom1010RecSFPS.hpp"
#include "ElectronRepulsionGeom1010RecSFSD.hpp"
#include "ElectronRepulsionGeom1010RecSFDS.hpp"
#include "ElectronRepulsionGeom1010RecSFSF.hpp"
#include "ElectronRepulsionGeom1010RecSFFS.hpp"
#include "ElectronRepulsionGeom1010RecSFPP.hpp"
#include "ElectronRepulsionGeom1010RecSFPD.hpp"
#include "ElectronRepulsionGeom1010RecSFDP.hpp"
#include "ElectronRepulsionGeom1010RecSFPF.hpp"
#include "ElectronRepulsionGeom1010RecSFFP.hpp"
#include "ElectronRepulsionGeom1010RecSFDD.hpp"
#include "ElectronRepulsionGeom1010RecSFDF.hpp"
#include "ElectronRepulsionGeom1010RecSFFD.hpp"
#include "ElectronRepulsionGeom1010RecSFFF.hpp"

#include "ElectronRepulsionGeom1010RecFSSS.hpp"
#include "ElectronRepulsionGeom1010RecFSSP.hpp"
#include "ElectronRepulsionGeom1010RecFSPS.hpp"
#include "ElectronRepulsionGeom1010RecFSSD.hpp"
#include "ElectronRepulsionGeom1010RecFSDS.hpp"
#include "ElectronRepulsionGeom1010RecFSSF.hpp"
#include "ElectronRepulsionGeom1010RecFSFS.hpp"
#include "ElectronRepulsionGeom1010RecFSPP.hpp"
#include "ElectronRepulsionGeom1010RecFSPD.hpp"
#include "ElectronRepulsionGeom1010RecFSDP.hpp"
#include "ElectronRepulsionGeom1010RecFSPF.hpp"
#include "ElectronRepulsionGeom1010RecFSFP.hpp"
#include "ElectronRepulsionGeom1010RecFSDD.hpp"
#include "ElectronRepulsionGeom1010RecFSDF.hpp"
#include "ElectronRepulsionGeom1010RecFSFD.hpp"
#include "ElectronRepulsionGeom1010RecFSFF.hpp"

#include "ElectronRepulsionGeom1010RecPPSS.hpp"
#include "ElectronRepulsionGeom1010RecPPSP.hpp"
#include "ElectronRepulsionGeom1010RecPPPS.hpp"
#include "ElectronRepulsionGeom1010RecPPSD.hpp"
#include "ElectronRepulsionGeom1010RecPPDS.hpp"
#include "ElectronRepulsionGeom1010RecPPSF.hpp"
#include "ElectronRepulsionGeom1010RecPPFS.hpp"
#include "ElectronRepulsionGeom1010RecPPPP.hpp"
#include "ElectronRepulsionGeom1010RecPPPD.hpp"
#include "ElectronRepulsionGeom1010RecPPDP.hpp"
#include "ElectronRepulsionGeom1010RecPPPF.hpp"
#include "ElectronRepulsionGeom1010RecPPFP.hpp"
#include "ElectronRepulsionGeom1010RecPPDD.hpp"
#include "ElectronRepulsionGeom1010RecPPDF.hpp"
#include "ElectronRepulsionGeom1010RecPPFD.hpp"
#include "ElectronRepulsionGeom1010RecPPFF.hpp"

#include "ElectronRepulsionGeom1010RecPDSS.hpp"
#include "ElectronRepulsionGeom1010RecPDSP.hpp"
#include "ElectronRepulsionGeom1010RecPDPS.hpp"
#include "ElectronRepulsionGeom1010RecPDSD.hpp"
#include "ElectronRepulsionGeom1010RecPDDS.hpp"
#include "ElectronRepulsionGeom1010RecPDSF.hpp"
#include "ElectronRepulsionGeom1010RecPDFS.hpp"
#include "ElectronRepulsionGeom1010RecPDPP.hpp"
#include "ElectronRepulsionGeom1010RecPDPD.hpp"
#include "ElectronRepulsionGeom1010RecPDDP.hpp"
#include "ElectronRepulsionGeom1010RecPDPF.hpp"
#include "ElectronRepulsionGeom1010RecPDFP.hpp"
#include "ElectronRepulsionGeom1010RecPDDD.hpp"
#include "ElectronRepulsionGeom1010RecPDDF.hpp"
#include "ElectronRepulsionGeom1010RecPDFD.hpp"
#include "ElectronRepulsionGeom1010RecPDFF.hpp"

#include "ElectronRepulsionGeom1010RecDPSS.hpp"
#include "ElectronRepulsionGeom1010RecDPSP.hpp"
#include "ElectronRepulsionGeom1010RecDPPS.hpp"
#include "ElectronRepulsionGeom1010RecDPSD.hpp"
#include "ElectronRepulsionGeom1010RecDPDS.hpp"
#include "ElectronRepulsionGeom1010RecDPSF.hpp"
#include "ElectronRepulsionGeom1010RecDPFS.hpp"
#include "ElectronRepulsionGeom1010RecDPPP.hpp"
#include "ElectronRepulsionGeom1010RecDPPD.hpp"
#include "ElectronRepulsionGeom1010RecDPDP.hpp"
#include "ElectronRepulsionGeom1010RecDPPF.hpp"
#include "ElectronRepulsionGeom1010RecDPFP.hpp"
#include "ElectronRepulsionGeom1010RecDPDD.hpp"
#include "ElectronRepulsionGeom1010RecDPDF.hpp"
#include "ElectronRepulsionGeom1010RecDPFD.hpp"
#include "ElectronRepulsionGeom1010RecDPFF.hpp"

#include "ElectronRepulsionGeom1010RecPFSS.hpp"
#include "ElectronRepulsionGeom1010RecPFSP.hpp"
#include "ElectronRepulsionGeom1010RecPFPS.hpp"
#include "ElectronRepulsionGeom1010RecPFSD.hpp"
#include "ElectronRepulsionGeom1010RecPFDS.hpp"
#include "ElectronRepulsionGeom1010RecPFSF.hpp"
#include "ElectronRepulsionGeom1010RecPFFS.hpp"
#include "ElectronRepulsionGeom1010RecPFPP.hpp"
#include "ElectronRepulsionGeom1010RecPFPD.hpp"
#include "ElectronRepulsionGeom1010RecPFDP.hpp"
#include "ElectronRepulsionGeom1010RecPFPF.hpp"
#include "ElectronRepulsionGeom1010RecPFFP.hpp"
#include "ElectronRepulsionGeom1010RecPFDD.hpp"
#include "ElectronRepulsionGeom1010RecPFDF.hpp"
#include "ElectronRepulsionGeom1010RecPFFD.hpp"
#include "ElectronRepulsionGeom1010RecPFFF.hpp"

#include "ElectronRepulsionGeom1010RecFPSS.hpp"
#include "ElectronRepulsionGeom1010RecFPSP.hpp"
#include "ElectronRepulsionGeom1010RecFPPS.hpp"
#include "ElectronRepulsionGeom1010RecFPSD.hpp"
#include "ElectronRepulsionGeom1010RecFPDS.hpp"
#include "ElectronRepulsionGeom1010RecFPSF.hpp"
#include "ElectronRepulsionGeom1010RecFPFS.hpp"
#include "ElectronRepulsionGeom1010RecFPPP.hpp"
#include "ElectronRepulsionGeom1010RecFPPD.hpp"
#include "ElectronRepulsionGeom1010RecFPDP.hpp"
#include "ElectronRepulsionGeom1010RecFPPF.hpp"
#include "ElectronRepulsionGeom1010RecFPFP.hpp"
#include "ElectronRepulsionGeom1010RecFPDD.hpp"
#include "ElectronRepulsionGeom1010RecFPDF.hpp"
#include "ElectronRepulsionGeom1010RecFPFD.hpp"
#include "ElectronRepulsionGeom1010RecFPFF.hpp"

#include "ElectronRepulsionGeom1010RecDDSS.hpp"
#include "ElectronRepulsionGeom1010RecDDSP.hpp"
#include "ElectronRepulsionGeom1010RecDDPS.hpp"
#include "ElectronRepulsionGeom1010RecDDSD.hpp"
#include "ElectronRepulsionGeom1010RecDDDS.hpp"
#include "ElectronRepulsionGeom1010RecDDSF.hpp"
#include "ElectronRepulsionGeom1010RecDDFS.hpp"
#include "ElectronRepulsionGeom1010RecDDPP.hpp"
#include "ElectronRepulsionGeom1010RecDDPD.hpp"
#include "ElectronRepulsionGeom1010RecDDDP.hpp"
#include "ElectronRepulsionGeom1010RecDDPF.hpp"
#include "ElectronRepulsionGeom1010RecDDFP.hpp"
#include "ElectronRepulsionGeom1010RecDDDD.hpp"
#include "ElectronRepulsionGeom1010RecDDDF.hpp"
#include "ElectronRepulsionGeom1010RecDDFD.hpp"
#include "ElectronRepulsionGeom1010RecDDFF.hpp"

#include "ElectronRepulsionGeom1010RecDFSS.hpp"
#include "ElectronRepulsionGeom1010RecDFSP.hpp"
#include "ElectronRepulsionGeom1010RecDFPS.hpp"
#include "ElectronRepulsionGeom1010RecDFSD.hpp"
#include "ElectronRepulsionGeom1010RecDFDS.hpp"
#include "ElectronRepulsionGeom1010RecDFSF.hpp"
#include "ElectronRepulsionGeom1010RecDFFS.hpp"
#include "ElectronRepulsionGeom1010RecDFPP.hpp"
#include "ElectronRepulsionGeom1010RecDFPD.hpp"
#include "ElectronRepulsionGeom1010RecDFDP.hpp"
#include "ElectronRepulsionGeom1010RecDFPF.hpp"
#include "ElectronRepulsionGeom1010RecDFFP.hpp"
#include "ElectronRepulsionGeom1010RecDFDD.hpp"
#include "ElectronRepulsionGeom1010RecDFDF.hpp"
#include "ElectronRepulsionGeom1010RecDFFD.hpp"
#include "ElectronRepulsionGeom1010RecDFFF.hpp"

#include "ElectronRepulsionGeom1010RecFDSS.hpp"
#include "ElectronRepulsionGeom1010RecFDSP.hpp"
#include "ElectronRepulsionGeom1010RecFDPS.hpp"
#include "ElectronRepulsionGeom1010RecFDSD.hpp"
#include "ElectronRepulsionGeom1010RecFDDS.hpp"
#include "ElectronRepulsionGeom1010RecFDSF.hpp"
#include "ElectronRepulsionGeom1010RecFDFS.hpp"
#include "ElectronRepulsionGeom1010RecFDPP.hpp"
#include "ElectronRepulsionGeom1010RecFDPD.hpp"
#include "ElectronRepulsionGeom1010RecFDDP.hpp"
#include "ElectronRepulsionGeom1010RecFDPF.hpp"
#include "ElectronRepulsionGeom1010RecFDFP.hpp"
#include "ElectronRepulsionGeom1010RecFDDD.hpp"
#include "ElectronRepulsionGeom1010RecFDDF.hpp"
#include "ElectronRepulsionGeom1010RecFDFD.hpp"
#include "ElectronRepulsionGeom1010RecFDFF.hpp"

#include "ElectronRepulsionGeom1010RecFFSS.hpp"
#include "ElectronRepulsionGeom1010RecFFSP.hpp"
#include "ElectronRepulsionGeom1010RecFFPS.hpp"
#include "ElectronRepulsionGeom1010RecFFSD.hpp"
#include "ElectronRepulsionGeom1010RecFFDS.hpp"
#include "ElectronRepulsionGeom1010RecFFSF.hpp"
#include "ElectronRepulsionGeom1010RecFFFS.hpp"
#include "ElectronRepulsionGeom1010RecFFPP.hpp"
#include "ElectronRepulsionGeom1010RecFFPD.hpp"
#include "ElectronRepulsionGeom1010RecFFDP.hpp"
#include "ElectronRepulsionGeom1010RecFFPF.hpp"
#include "ElectronRepulsionGeom1010RecFFFP.hpp"
#include "ElectronRepulsionGeom1010RecFFDD.hpp"
#include "ElectronRepulsionGeom1010RecFFDF.hpp"
#include "ElectronRepulsionGeom1010RecFFFD.hpp"
#include "ElectronRepulsionGeom1010RecFFFF.hpp"

#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
compute_geom_1010(T&                               distributor,
                  const CGtoPairBlock&             bra_gto_pair_block,
                  const CGtoPairBlock&             ket_gto_pair_block,
                  const std::pair<size_t, size_t>& bra_indices,
                  const std::pair<size_t, size_t>& ket_indices) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();
    
    // leading [SS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
   
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ssdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ssfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ssfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_spdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_spfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_spfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_psdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_psfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_psdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_psfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_psff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sdfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sdfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dsdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dsfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dsfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sfps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sfds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sffs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sfdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sffp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sffd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fsps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fsds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fsfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fsdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fsfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fsfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ppsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ppdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ppfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ppdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ppfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ppff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pdfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pdfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dpdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dpfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dpfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pfps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pfds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pffs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pfdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pffp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pffd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_pfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fpps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fpds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fpfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fpdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fpfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fpfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ddsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ddpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ddfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ddfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ddff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dfps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dfds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dffs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dfdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dffp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dffd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_dfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fdps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fdds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fdfs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fdfp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fdfd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_fdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ffss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ffsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ffps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ffsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ffds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ffsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_fffs(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ffpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ffpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ffdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ffpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_fffp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ffdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ffdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_fffd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1010_ffff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;

    std::abort();
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1100Func_hpp */
