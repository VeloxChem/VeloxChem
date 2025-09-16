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

#ifndef ElectronRepulsionGeom1100Func_hpp
#define ElectronRepulsionGeom1100Func_hpp

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom1100RecSSSS.hpp"
#include "ElectronRepulsionGeom1100RecSSSP.hpp"
#include "ElectronRepulsionGeom1100RecSSSD.hpp"
#include "ElectronRepulsionGeom1100RecSSSF.hpp"
#include "ElectronRepulsionGeom1100RecSSPP.hpp"
#include "ElectronRepulsionGeom1100RecSSPD.hpp"
#include "ElectronRepulsionGeom1100RecSSPF.hpp"
#include "ElectronRepulsionGeom1100RecSSDD.hpp"
#include "ElectronRepulsionGeom1100RecSSDF.hpp"
#include "ElectronRepulsionGeom1100RecSSFF.hpp"
#include "ElectronRepulsionGeom1100RecSPSS.hpp"
#include "ElectronRepulsionGeom1100RecSPSP.hpp"
#include "ElectronRepulsionGeom1100RecSPSD.hpp"
#include "ElectronRepulsionGeom1100RecSPSF.hpp"
#include "ElectronRepulsionGeom1100RecSPPP.hpp"
#include "ElectronRepulsionGeom1100RecSPPD.hpp"
#include "ElectronRepulsionGeom1100RecSPPF.hpp"
#include "ElectronRepulsionGeom1100RecSPDD.hpp"
#include "ElectronRepulsionGeom1100RecSPDF.hpp"
#include "ElectronRepulsionGeom1100RecSPFF.hpp"
#include "ElectronRepulsionGeom1100RecSDSS.hpp"
#include "ElectronRepulsionGeom1100RecSDSP.hpp"
#include "ElectronRepulsionGeom1100RecSDSD.hpp"
#include "ElectronRepulsionGeom1100RecSDSF.hpp"
#include "ElectronRepulsionGeom1100RecSDPP.hpp"
#include "ElectronRepulsionGeom1100RecSDPD.hpp"
#include "ElectronRepulsionGeom1100RecSDPF.hpp"
#include "ElectronRepulsionGeom1100RecSDDD.hpp"
#include "ElectronRepulsionGeom1100RecSDDF.hpp"
#include "ElectronRepulsionGeom1100RecSDFF.hpp"
#include "ElectronRepulsionGeom1100RecSFSS.hpp"
#include "ElectronRepulsionGeom1100RecSFSP.hpp"
#include "ElectronRepulsionGeom1100RecSFSD.hpp"
#include "ElectronRepulsionGeom1100RecSFSF.hpp"
#include "ElectronRepulsionGeom1100RecSFPP.hpp"
#include "ElectronRepulsionGeom1100RecSFPD.hpp"
#include "ElectronRepulsionGeom1100RecSFPF.hpp"
#include "ElectronRepulsionGeom1100RecSFDD.hpp"
#include "ElectronRepulsionGeom1100RecSFDF.hpp"
#include "ElectronRepulsionGeom1100RecSFFF.hpp"
#include "ElectronRepulsionGeom1100RecPSSS.hpp"
#include "ElectronRepulsionGeom1100RecPSSP.hpp"
#include "ElectronRepulsionGeom1100RecPSSD.hpp"
#include "ElectronRepulsionGeom1100RecPSSF.hpp"
#include "ElectronRepulsionGeom1100RecPSPP.hpp"
#include "ElectronRepulsionGeom1100RecPSPD.hpp"
#include "ElectronRepulsionGeom1100RecPSPF.hpp"
#include "ElectronRepulsionGeom1100RecPSDD.hpp"
#include "ElectronRepulsionGeom1100RecPSDF.hpp"
#include "ElectronRepulsionGeom1100RecPSFF.hpp"
#include "ElectronRepulsionGeom1100RecDSSS.hpp"
#include "ElectronRepulsionGeom1100RecDSSP.hpp"
#include "ElectronRepulsionGeom1100RecDSSD.hpp"
#include "ElectronRepulsionGeom1100RecDSSF.hpp"
#include "ElectronRepulsionGeom1100RecDSPP.hpp"
#include "ElectronRepulsionGeom1100RecDSPD.hpp"
#include "ElectronRepulsionGeom1100RecDSPF.hpp"
#include "ElectronRepulsionGeom1100RecDSDD.hpp"
#include "ElectronRepulsionGeom1100RecDSDF.hpp"
#include "ElectronRepulsionGeom1100RecDSFF.hpp"
#include "ElectronRepulsionGeom1100RecFSSS.hpp"
#include "ElectronRepulsionGeom1100RecFSSP.hpp"
#include "ElectronRepulsionGeom1100RecFSSD.hpp"
#include "ElectronRepulsionGeom1100RecFSSF.hpp"
#include "ElectronRepulsionGeom1100RecFSPP.hpp"
#include "ElectronRepulsionGeom1100RecFSPD.hpp"
#include "ElectronRepulsionGeom1100RecFSPF.hpp"
#include "ElectronRepulsionGeom1100RecFSDD.hpp"
#include "ElectronRepulsionGeom1100RecFSDF.hpp"
#include "ElectronRepulsionGeom1100RecFSFF.hpp"
#include "ElectronRepulsionGeom1100RecPPSS.hpp"
#include "ElectronRepulsionGeom1100RecPPSP.hpp"
#include "ElectronRepulsionGeom1100RecPPSD.hpp"
#include "ElectronRepulsionGeom1100RecPPSF.hpp"
#include "ElectronRepulsionGeom1100RecPPPP.hpp"
#include "ElectronRepulsionGeom1100RecPPPD.hpp"
#include "ElectronRepulsionGeom1100RecPPPF.hpp"
#include "ElectronRepulsionGeom1100RecPPDD.hpp"
#include "ElectronRepulsionGeom1100RecPPDF.hpp"
#include "ElectronRepulsionGeom1100RecPPFF.hpp"
#include "ElectronRepulsionGeom1100RecPDSS.hpp"
#include "ElectronRepulsionGeom1100RecPDSP.hpp"
#include "ElectronRepulsionGeom1100RecPDSD.hpp"
#include "ElectronRepulsionGeom1100RecPDSF.hpp"
#include "ElectronRepulsionGeom1100RecPDPP.hpp"
#include "ElectronRepulsionGeom1100RecPDPD.hpp"
#include "ElectronRepulsionGeom1100RecPDPF.hpp"
#include "ElectronRepulsionGeom1100RecPDDD.hpp"
#include "ElectronRepulsionGeom1100RecPDDF.hpp"
#include "ElectronRepulsionGeom1100RecPDFF.hpp"
#include "ElectronRepulsionGeom1100RecPFSS.hpp"
#include "ElectronRepulsionGeom1100RecPFSP.hpp"
#include "ElectronRepulsionGeom1100RecPFSD.hpp"
#include "ElectronRepulsionGeom1100RecPFSF.hpp"
#include "ElectronRepulsionGeom1100RecPFPP.hpp"
#include "ElectronRepulsionGeom1100RecPFPD.hpp"
#include "ElectronRepulsionGeom1100RecPFPF.hpp"
#include "ElectronRepulsionGeom1100RecPFDD.hpp"
#include "ElectronRepulsionGeom1100RecPFDF.hpp"
#include "ElectronRepulsionGeom1100RecPFFF.hpp"
#include "ElectronRepulsionGeom1100RecDPSS.hpp"
#include "ElectronRepulsionGeom1100RecDPSP.hpp"
#include "ElectronRepulsionGeom1100RecDPSD.hpp"
#include "ElectronRepulsionGeom1100RecDPSF.hpp"
#include "ElectronRepulsionGeom1100RecDPPP.hpp"
#include "ElectronRepulsionGeom1100RecDPPD.hpp"
#include "ElectronRepulsionGeom1100RecDPPF.hpp"
#include "ElectronRepulsionGeom1100RecDPDD.hpp"
#include "ElectronRepulsionGeom1100RecDPDF.hpp"
#include "ElectronRepulsionGeom1100RecDPFF.hpp"
#include "ElectronRepulsionGeom1100RecFPSS.hpp"
#include "ElectronRepulsionGeom1100RecFPSP.hpp"
#include "ElectronRepulsionGeom1100RecFPSD.hpp"
#include "ElectronRepulsionGeom1100RecFPSF.hpp"
#include "ElectronRepulsionGeom1100RecFPPP.hpp"
#include "ElectronRepulsionGeom1100RecFPPD.hpp"
#include "ElectronRepulsionGeom1100RecFPPF.hpp"
#include "ElectronRepulsionGeom1100RecFPDD.hpp"
#include "ElectronRepulsionGeom1100RecFPDF.hpp"
#include "ElectronRepulsionGeom1100RecFPFF.hpp"
#include "ElectronRepulsionGeom1100RecDDSS.hpp"
#include "ElectronRepulsionGeom1100RecDDSP.hpp"
#include "ElectronRepulsionGeom1100RecDDSD.hpp"
#include "ElectronRepulsionGeom1100RecDDSF.hpp"
#include "ElectronRepulsionGeom1100RecDDPP.hpp"
#include "ElectronRepulsionGeom1100RecDDPD.hpp"
#include "ElectronRepulsionGeom1100RecDDPF.hpp"
#include "ElectronRepulsionGeom1100RecDDDD.hpp"
#include "ElectronRepulsionGeom1100RecDDDF.hpp"
#include "ElectronRepulsionGeom1100RecDDFF.hpp"
#include "ElectronRepulsionGeom1100RecDFSS.hpp"
#include "ElectronRepulsionGeom1100RecDFSP.hpp"
#include "ElectronRepulsionGeom1100RecDFSD.hpp"
#include "ElectronRepulsionGeom1100RecDFSF.hpp"
#include "ElectronRepulsionGeom1100RecDFPP.hpp"
#include "ElectronRepulsionGeom1100RecDFPD.hpp"
#include "ElectronRepulsionGeom1100RecDFPF.hpp"
#include "ElectronRepulsionGeom1100RecDFDD.hpp"
#include "ElectronRepulsionGeom1100RecDFDF.hpp"
#include "ElectronRepulsionGeom1100RecDFFF.hpp"
#include "ElectronRepulsionGeom1100RecFDSS.hpp"
#include "ElectronRepulsionGeom1100RecFDSP.hpp"
#include "ElectronRepulsionGeom1100RecFDSD.hpp"
#include "ElectronRepulsionGeom1100RecFDSF.hpp"
#include "ElectronRepulsionGeom1100RecFDPP.hpp"
#include "ElectronRepulsionGeom1100RecFDPD.hpp"
#include "ElectronRepulsionGeom1100RecFDPF.hpp"
#include "ElectronRepulsionGeom1100RecFDDD.hpp"
#include "ElectronRepulsionGeom1100RecFDDF.hpp"
#include "ElectronRepulsionGeom1100RecFDFF.hpp"
#include "ElectronRepulsionGeom1100RecFFSS.hpp"
#include "ElectronRepulsionGeom1100RecFFSP.hpp"
#include "ElectronRepulsionGeom1100RecFFSD.hpp"
#include "ElectronRepulsionGeom1100RecFFSF.hpp"
#include "ElectronRepulsionGeom1100RecFFPP.hpp"
#include "ElectronRepulsionGeom1100RecFFPD.hpp"
#include "ElectronRepulsionGeom1100RecFFPF.hpp"
#include "ElectronRepulsionGeom1100RecFFDD.hpp"
#include "ElectronRepulsionGeom1100RecFFDF.hpp"
#include "ElectronRepulsionGeom1100RecFFFF.hpp"

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
compute_geom_1100(T&                               distributor,
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
        erirec::comp_electron_repulsion_geom1100_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_psdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_psff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_dsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_fsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ppsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ppdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ppff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_pfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_pfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_pfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_pfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_dpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_fpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ddsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ddpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ddff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_dfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_dfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_dfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_dfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_fdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_fdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_fddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_fdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1100_ffss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_ffsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ffsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ffsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1100_ffpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ffpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ffpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1100_ffdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ffdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1100_ffff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;

    std::abort();
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1100Func_hpp */
