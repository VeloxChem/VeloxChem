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

#ifndef ElectronRepulsionGeom2000Func_hpp
#define ElectronRepulsionGeom2000Func_hpp

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom2000RecSSSS.hpp"
#include "ElectronRepulsionGeom2000RecSSSP.hpp"
#include "ElectronRepulsionGeom2000RecSSSD.hpp"
#include "ElectronRepulsionGeom2000RecSSSF.hpp"
#include "ElectronRepulsionGeom2000RecSSPP.hpp"
#include "ElectronRepulsionGeom2000RecSSPD.hpp"
#include "ElectronRepulsionGeom2000RecSSPF.hpp"
#include "ElectronRepulsionGeom2000RecSSDD.hpp"
#include "ElectronRepulsionGeom2000RecSSDF.hpp"
#include "ElectronRepulsionGeom2000RecSSFF.hpp"
#include "ElectronRepulsionGeom2000RecSPSS.hpp"
#include "ElectronRepulsionGeom2000RecSPSP.hpp"
#include "ElectronRepulsionGeom2000RecSPSD.hpp"
#include "ElectronRepulsionGeom2000RecSPSF.hpp"
#include "ElectronRepulsionGeom2000RecSPPP.hpp"
#include "ElectronRepulsionGeom2000RecSPPD.hpp"
#include "ElectronRepulsionGeom2000RecSPPF.hpp"
#include "ElectronRepulsionGeom2000RecSPDD.hpp"
#include "ElectronRepulsionGeom2000RecSPDF.hpp"
#include "ElectronRepulsionGeom2000RecSPFF.hpp"
#include "ElectronRepulsionGeom2000RecSDSS.hpp"
#include "ElectronRepulsionGeom2000RecSDSP.hpp"
#include "ElectronRepulsionGeom2000RecSDSD.hpp"
#include "ElectronRepulsionGeom2000RecSDSF.hpp"
#include "ElectronRepulsionGeom2000RecSDPP.hpp"
#include "ElectronRepulsionGeom2000RecSDPD.hpp"
#include "ElectronRepulsionGeom2000RecSDPF.hpp"
#include "ElectronRepulsionGeom2000RecSDDD.hpp"
#include "ElectronRepulsionGeom2000RecSDDF.hpp"
#include "ElectronRepulsionGeom2000RecSDFF.hpp"
#include "ElectronRepulsionGeom2000RecSFSS.hpp"
#include "ElectronRepulsionGeom2000RecSFSP.hpp"
#include "ElectronRepulsionGeom2000RecSFSD.hpp"
#include "ElectronRepulsionGeom2000RecSFSF.hpp"
#include "ElectronRepulsionGeom2000RecSFPP.hpp"
#include "ElectronRepulsionGeom2000RecSFPD.hpp"
#include "ElectronRepulsionGeom2000RecSFPF.hpp"
#include "ElectronRepulsionGeom2000RecSFDD.hpp"
#include "ElectronRepulsionGeom2000RecSFDF.hpp"
#include "ElectronRepulsionGeom2000RecSFFF.hpp"
#include "ElectronRepulsionGeom2000RecPSSS.hpp"
#include "ElectronRepulsionGeom2000RecPSSP.hpp"
#include "ElectronRepulsionGeom2000RecPSSD.hpp"
#include "ElectronRepulsionGeom2000RecPSSF.hpp"
#include "ElectronRepulsionGeom2000RecPSPP.hpp"
#include "ElectronRepulsionGeom2000RecPSPD.hpp"
#include "ElectronRepulsionGeom2000RecPSPF.hpp"
#include "ElectronRepulsionGeom2000RecPSDD.hpp"
#include "ElectronRepulsionGeom2000RecPSDF.hpp"
#include "ElectronRepulsionGeom2000RecPSFF.hpp"
#include "ElectronRepulsionGeom2000RecDSSS.hpp"
#include "ElectronRepulsionGeom2000RecDSSP.hpp"
#include "ElectronRepulsionGeom2000RecDSSD.hpp"
#include "ElectronRepulsionGeom2000RecDSSF.hpp"
#include "ElectronRepulsionGeom2000RecDSPP.hpp"
#include "ElectronRepulsionGeom2000RecDSPD.hpp"
#include "ElectronRepulsionGeom2000RecDSPF.hpp"
#include "ElectronRepulsionGeom2000RecDSDD.hpp"
#include "ElectronRepulsionGeom2000RecDSDF.hpp"
#include "ElectronRepulsionGeom2000RecDSFF.hpp"
#include "ElectronRepulsionGeom2000RecFSSS.hpp"
#include "ElectronRepulsionGeom2000RecFSSP.hpp"
#include "ElectronRepulsionGeom2000RecFSSD.hpp"
#include "ElectronRepulsionGeom2000RecFSSF.hpp"
#include "ElectronRepulsionGeom2000RecFSPP.hpp"
#include "ElectronRepulsionGeom2000RecFSPD.hpp"
#include "ElectronRepulsionGeom2000RecFSPF.hpp"
#include "ElectronRepulsionGeom2000RecFSDD.hpp"
#include "ElectronRepulsionGeom2000RecFSDF.hpp"
#include "ElectronRepulsionGeom2000RecFSFF.hpp"
#include "ElectronRepulsionGeom2000RecPPSS.hpp"
#include "ElectronRepulsionGeom2000RecPPSP.hpp"
#include "ElectronRepulsionGeom2000RecPPSD.hpp"
#include "ElectronRepulsionGeom2000RecPPSF.hpp"
#include "ElectronRepulsionGeom2000RecPPPP.hpp"
#include "ElectronRepulsionGeom2000RecPPPD.hpp"
#include "ElectronRepulsionGeom2000RecPPPF.hpp"
#include "ElectronRepulsionGeom2000RecPPDD.hpp"
#include "ElectronRepulsionGeom2000RecPPDF.hpp"
#include "ElectronRepulsionGeom2000RecPPFF.hpp"
#include "ElectronRepulsionGeom2000RecPDSS.hpp"
#include "ElectronRepulsionGeom2000RecPDSP.hpp"
#include "ElectronRepulsionGeom2000RecPDSD.hpp"
#include "ElectronRepulsionGeom2000RecPDSF.hpp"
#include "ElectronRepulsionGeom2000RecPDPP.hpp"
#include "ElectronRepulsionGeom2000RecPDPD.hpp"
#include "ElectronRepulsionGeom2000RecPDPF.hpp"
#include "ElectronRepulsionGeom2000RecPDDD.hpp"
#include "ElectronRepulsionGeom2000RecPDDF.hpp"
#include "ElectronRepulsionGeom2000RecPDFF.hpp"
#include "ElectronRepulsionGeom2000RecPFSS.hpp"
#include "ElectronRepulsionGeom2000RecPFSP.hpp"
#include "ElectronRepulsionGeom2000RecPFSD.hpp"
#include "ElectronRepulsionGeom2000RecPFSF.hpp"
#include "ElectronRepulsionGeom2000RecPFPP.hpp"
#include "ElectronRepulsionGeom2000RecPFPD.hpp"
#include "ElectronRepulsionGeom2000RecPFPF.hpp"
#include "ElectronRepulsionGeom2000RecPFDD.hpp"
#include "ElectronRepulsionGeom2000RecPFDF.hpp"
#include "ElectronRepulsionGeom2000RecPFFF.hpp"
#include "ElectronRepulsionGeom2000RecDPSS.hpp"
#include "ElectronRepulsionGeom2000RecDPSP.hpp"
#include "ElectronRepulsionGeom2000RecDPSD.hpp"
#include "ElectronRepulsionGeom2000RecDPSF.hpp"
#include "ElectronRepulsionGeom2000RecDPPP.hpp"
#include "ElectronRepulsionGeom2000RecDPPD.hpp"
#include "ElectronRepulsionGeom2000RecDPPF.hpp"
#include "ElectronRepulsionGeom2000RecDPDD.hpp"
#include "ElectronRepulsionGeom2000RecDPDF.hpp"
#include "ElectronRepulsionGeom2000RecDPFF.hpp"
#include "ElectronRepulsionGeom2000RecFPSS.hpp"
#include "ElectronRepulsionGeom2000RecFPSP.hpp"
#include "ElectronRepulsionGeom2000RecFPSD.hpp"
#include "ElectronRepulsionGeom2000RecFPSF.hpp"
#include "ElectronRepulsionGeom2000RecFPPP.hpp"
#include "ElectronRepulsionGeom2000RecFPPD.hpp"
#include "ElectronRepulsionGeom2000RecFPPF.hpp"
#include "ElectronRepulsionGeom2000RecFPDD.hpp"
#include "ElectronRepulsionGeom2000RecFPDF.hpp"
#include "ElectronRepulsionGeom2000RecFPFF.hpp"
#include "ElectronRepulsionGeom2000RecDDSS.hpp"
#include "ElectronRepulsionGeom2000RecDDSP.hpp"
#include "ElectronRepulsionGeom2000RecDDSD.hpp"
#include "ElectronRepulsionGeom2000RecDDSF.hpp"
#include "ElectronRepulsionGeom2000RecDDPP.hpp"
#include "ElectronRepulsionGeom2000RecDDPD.hpp"
#include "ElectronRepulsionGeom2000RecDDPF.hpp"
#include "ElectronRepulsionGeom2000RecDDDD.hpp"
#include "ElectronRepulsionGeom2000RecDDDF.hpp"
#include "ElectronRepulsionGeom2000RecDDFF.hpp"
#include "ElectronRepulsionGeom2000RecDFSS.hpp"
#include "ElectronRepulsionGeom2000RecDFSP.hpp"
#include "ElectronRepulsionGeom2000RecDFSD.hpp"
#include "ElectronRepulsionGeom2000RecDFSF.hpp"
#include "ElectronRepulsionGeom2000RecDFPP.hpp"
#include "ElectronRepulsionGeom2000RecDFPD.hpp"
#include "ElectronRepulsionGeom2000RecDFPF.hpp"
#include "ElectronRepulsionGeom2000RecDFDD.hpp"
#include "ElectronRepulsionGeom2000RecDFDF.hpp"
#include "ElectronRepulsionGeom2000RecDFFF.hpp"
#include "ElectronRepulsionGeom2000RecFDSS.hpp"
#include "ElectronRepulsionGeom2000RecFDSP.hpp"
#include "ElectronRepulsionGeom2000RecFDSD.hpp"
#include "ElectronRepulsionGeom2000RecFDSF.hpp"
#include "ElectronRepulsionGeom2000RecFDPP.hpp"
#include "ElectronRepulsionGeom2000RecFDPD.hpp"
#include "ElectronRepulsionGeom2000RecFDPF.hpp"
#include "ElectronRepulsionGeom2000RecFDDD.hpp"
#include "ElectronRepulsionGeom2000RecFDDF.hpp"
#include "ElectronRepulsionGeom2000RecFDFF.hpp"
#include "ElectronRepulsionGeom2000RecFFSS.hpp"
#include "ElectronRepulsionGeom2000RecFFSP.hpp"
#include "ElectronRepulsionGeom2000RecFFSD.hpp"
#include "ElectronRepulsionGeom2000RecFFSF.hpp"
#include "ElectronRepulsionGeom2000RecFFPP.hpp"
#include "ElectronRepulsionGeom2000RecFFPD.hpp"
#include "ElectronRepulsionGeom2000RecFFPF.hpp"
#include "ElectronRepulsionGeom2000RecFFDD.hpp"
#include "ElectronRepulsionGeom2000RecFFDF.hpp"
#include "ElectronRepulsionGeom2000RecFFFF.hpp"

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
compute_geom_2000(T&                               distributor,
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
        erirec::comp_electron_repulsion_geom2000_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_psdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_psff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_dsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_fsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ppsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ppdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ppff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_pfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_pfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_pfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_pfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_dpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_fpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ddsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ddpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ddff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_dfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_dfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_dfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_dfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_fdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_fdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_fddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_fdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_ffss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_ffsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ffsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ffsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom2000_ffpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ffpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ffpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom2000_ffdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ffdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom2000_ffff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;

    std::abort();
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom2000Func_hpp */
