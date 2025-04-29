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

#ifndef ElectronRepulsionGeom1000Func_hpp
#define ElectronRepulsionGeom1000Func_hpp

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom1000RecSSSS.hpp"
#include "ElectronRepulsionGeom1000RecSSSP.hpp"
#include "ElectronRepulsionGeom1000RecSSSD.hpp"
#include "ElectronRepulsionGeom1000RecSSSF.hpp"
#include "ElectronRepulsionGeom1000RecSSPP.hpp"
#include "ElectronRepulsionGeom1000RecSSPD.hpp"
#include "ElectronRepulsionGeom1000RecSSPF.hpp"
#include "ElectronRepulsionGeom1000RecSSDD.hpp"
#include "ElectronRepulsionGeom1000RecSSDF.hpp"
#include "ElectronRepulsionGeom1000RecSSFF.hpp"
#include "ElectronRepulsionGeom1000RecSPSS.hpp"
#include "ElectronRepulsionGeom1000RecSPSP.hpp"
#include "ElectronRepulsionGeom1000RecSPSD.hpp"
#include "ElectronRepulsionGeom1000RecSPSF.hpp"
#include "ElectronRepulsionGeom1000RecSPPP.hpp"
#include "ElectronRepulsionGeom1000RecSPPD.hpp"
#include "ElectronRepulsionGeom1000RecSPPF.hpp"
#include "ElectronRepulsionGeom1000RecSPDD.hpp"
#include "ElectronRepulsionGeom1000RecSPDF.hpp"
#include "ElectronRepulsionGeom1000RecSPFF.hpp"
#include "ElectronRepulsionGeom1000RecPSSS.hpp"
#include "ElectronRepulsionGeom1000RecPSSP.hpp"
#include "ElectronRepulsionGeom1000RecPSSD.hpp"
#include "ElectronRepulsionGeom1000RecPSSF.hpp"
#include "ElectronRepulsionGeom1000RecPSPP.hpp"
#include "ElectronRepulsionGeom1000RecPSPD.hpp"
#include "ElectronRepulsionGeom1000RecPSPF.hpp"
#include "ElectronRepulsionGeom1000RecPSDD.hpp"
#include "ElectronRepulsionGeom1000RecPSDF.hpp"
#include "ElectronRepulsionGeom1000RecPSFF.hpp"
#include "ElectronRepulsionGeom1000RecPPSS.hpp"
#include "ElectronRepulsionGeom1000RecPPSP.hpp"
#include "ElectronRepulsionGeom1000RecPPSD.hpp"
#include "ElectronRepulsionGeom1000RecPPSF.hpp"
#include "ElectronRepulsionGeom1000RecPPPP.hpp"
#include "ElectronRepulsionGeom1000RecPPPD.hpp"
#include "ElectronRepulsionGeom1000RecPPPF.hpp"
#include "ElectronRepulsionGeom1000RecPPDD.hpp"
#include "ElectronRepulsionGeom1000RecPPDF.hpp"
#include "ElectronRepulsionGeom1000RecPPFF.hpp"
#include "ElectronRepulsionGeom1000RecSDSS.hpp"
#include "ElectronRepulsionGeom1000RecSDSP.hpp"
#include "ElectronRepulsionGeom1000RecSDSD.hpp"
#include "ElectronRepulsionGeom1000RecSDSF.hpp"
#include "ElectronRepulsionGeom1000RecSDPP.hpp"
#include "ElectronRepulsionGeom1000RecSDPD.hpp"
#include "ElectronRepulsionGeom1000RecSDPF.hpp"
#include "ElectronRepulsionGeom1000RecSDDD.hpp"
#include "ElectronRepulsionGeom1000RecSDDF.hpp"
#include "ElectronRepulsionGeom1000RecSDFF.hpp"
#include "ElectronRepulsionGeom1000RecDSSS.hpp"
#include "ElectronRepulsionGeom1000RecDSSP.hpp"
#include "ElectronRepulsionGeom1000RecDSSD.hpp"
#include "ElectronRepulsionGeom1000RecDSSF.hpp"
#include "ElectronRepulsionGeom1000RecDSPP.hpp"
#include "ElectronRepulsionGeom1000RecDSPD.hpp"
#include "ElectronRepulsionGeom1000RecDSPF.hpp"
#include "ElectronRepulsionGeom1000RecDSDD.hpp"
#include "ElectronRepulsionGeom1000RecDSDF.hpp"
#include "ElectronRepulsionGeom1000RecDSFF.hpp"
#include "ElectronRepulsionGeom1000RecPDSS.hpp"
#include "ElectronRepulsionGeom1000RecPDSP.hpp"
#include "ElectronRepulsionGeom1000RecPDSD.hpp"
#include "ElectronRepulsionGeom1000RecPDSF.hpp"
#include "ElectronRepulsionGeom1000RecPDPP.hpp"
#include "ElectronRepulsionGeom1000RecPDPD.hpp"
#include "ElectronRepulsionGeom1000RecPDPF.hpp"
#include "ElectronRepulsionGeom1000RecPDDD.hpp"
#include "ElectronRepulsionGeom1000RecPDDF.hpp"
#include "ElectronRepulsionGeom1000RecPDFF.hpp"
#include "ElectronRepulsionGeom1000RecDPSS.hpp"
#include "ElectronRepulsionGeom1000RecDPSP.hpp"
#include "ElectronRepulsionGeom1000RecDPSD.hpp"
#include "ElectronRepulsionGeom1000RecDPSF.hpp"
#include "ElectronRepulsionGeom1000RecDPPP.hpp"
#include "ElectronRepulsionGeom1000RecDPPD.hpp"
#include "ElectronRepulsionGeom1000RecDPPF.hpp"
#include "ElectronRepulsionGeom1000RecDPDD.hpp"
#include "ElectronRepulsionGeom1000RecDPDF.hpp"
#include "ElectronRepulsionGeom1000RecDPFF.hpp"
#include "ElectronRepulsionGeom1000RecSFSS.hpp"
#include "ElectronRepulsionGeom1000RecSFSP.hpp"
#include "ElectronRepulsionGeom1000RecSFSD.hpp"
#include "ElectronRepulsionGeom1000RecSFSF.hpp"
#include "ElectronRepulsionGeom1000RecSFPP.hpp"
#include "ElectronRepulsionGeom1000RecSFPD.hpp"
#include "ElectronRepulsionGeom1000RecSFPF.hpp"
#include "ElectronRepulsionGeom1000RecSFDD.hpp"
#include "ElectronRepulsionGeom1000RecSFDF.hpp"
#include "ElectronRepulsionGeom1000RecSFFF.hpp"
#include "ElectronRepulsionGeom1000RecFSSS.hpp"
#include "ElectronRepulsionGeom1000RecFSSP.hpp"
#include "ElectronRepulsionGeom1000RecFSSD.hpp"
#include "ElectronRepulsionGeom1000RecFSSF.hpp"
#include "ElectronRepulsionGeom1000RecFSPP.hpp"
#include "ElectronRepulsionGeom1000RecFSPD.hpp"
#include "ElectronRepulsionGeom1000RecFSPF.hpp"
#include "ElectronRepulsionGeom1000RecFSDD.hpp"
#include "ElectronRepulsionGeom1000RecFSDF.hpp"
#include "ElectronRepulsionGeom1000RecFSFF.hpp"
#include "ElectronRepulsionGeom1000RecDDSS.hpp"
#include "ElectronRepulsionGeom1000RecDDSP.hpp"
#include "ElectronRepulsionGeom1000RecDDSD.hpp"
#include "ElectronRepulsionGeom1000RecDDSF.hpp"
#include "ElectronRepulsionGeom1000RecDDPP.hpp"
#include "ElectronRepulsionGeom1000RecDDPD.hpp"
#include "ElectronRepulsionGeom1000RecDDPF.hpp"
#include "ElectronRepulsionGeom1000RecDDDD.hpp"
#include "ElectronRepulsionGeom1000RecDDDF.hpp"
#include "ElectronRepulsionGeom1000RecDDFF.hpp"
#include "ElectronRepulsionGeom1000RecPFSS.hpp"
#include "ElectronRepulsionGeom1000RecPFSP.hpp"
#include "ElectronRepulsionGeom1000RecPFSD.hpp"
#include "ElectronRepulsionGeom1000RecPFSF.hpp"
#include "ElectronRepulsionGeom1000RecPFPP.hpp"
#include "ElectronRepulsionGeom1000RecPFPD.hpp"
#include "ElectronRepulsionGeom1000RecPFPF.hpp"
#include "ElectronRepulsionGeom1000RecPFDD.hpp"
#include "ElectronRepulsionGeom1000RecPFDF.hpp"
#include "ElectronRepulsionGeom1000RecPFFF.hpp"
#include "ElectronRepulsionGeom1000RecFPSS.hpp"
#include "ElectronRepulsionGeom1000RecFPSP.hpp"
#include "ElectronRepulsionGeom1000RecFPSD.hpp"
#include "ElectronRepulsionGeom1000RecFPSF.hpp"
#include "ElectronRepulsionGeom1000RecFPPP.hpp"
#include "ElectronRepulsionGeom1000RecFPPD.hpp"
#include "ElectronRepulsionGeom1000RecFPPF.hpp"
#include "ElectronRepulsionGeom1000RecFPDD.hpp"
#include "ElectronRepulsionGeom1000RecFPDF.hpp"
#include "ElectronRepulsionGeom1000RecFPFF.hpp"
#include "ElectronRepulsionGeom1000RecDFSS.hpp"
#include "ElectronRepulsionGeom1000RecDFSP.hpp"
#include "ElectronRepulsionGeom1000RecDFSD.hpp"
#include "ElectronRepulsionGeom1000RecDFSF.hpp"
#include "ElectronRepulsionGeom1000RecDFPP.hpp"
#include "ElectronRepulsionGeom1000RecDFPD.hpp"
#include "ElectronRepulsionGeom1000RecDFPF.hpp"
#include "ElectronRepulsionGeom1000RecDFDD.hpp"
#include "ElectronRepulsionGeom1000RecDFDF.hpp"
#include "ElectronRepulsionGeom1000RecDFFF.hpp"
#include "ElectronRepulsionGeom1000RecFDSS.hpp"
#include "ElectronRepulsionGeom1000RecFDSP.hpp"
#include "ElectronRepulsionGeom1000RecFDSD.hpp"
#include "ElectronRepulsionGeom1000RecFDSF.hpp"
#include "ElectronRepulsionGeom1000RecFDPP.hpp"
#include "ElectronRepulsionGeom1000RecFDPD.hpp"
#include "ElectronRepulsionGeom1000RecFDPF.hpp"
#include "ElectronRepulsionGeom1000RecFDDD.hpp"
#include "ElectronRepulsionGeom1000RecFDDF.hpp"
#include "ElectronRepulsionGeom1000RecFDFF.hpp"
#include "ElectronRepulsionGeom1000RecFFSS.hpp"
#include "ElectronRepulsionGeom1000RecFFSP.hpp"
#include "ElectronRepulsionGeom1000RecFFSD.hpp"
#include "ElectronRepulsionGeom1000RecFFSF.hpp"
#include "ElectronRepulsionGeom1000RecFFPP.hpp"
#include "ElectronRepulsionGeom1000RecFFPD.hpp"
#include "ElectronRepulsionGeom1000RecFFPF.hpp"
#include "ElectronRepulsionGeom1000RecFFDD.hpp"
#include "ElectronRepulsionGeom1000RecFFDF.hpp"
#include "ElectronRepulsionGeom1000RecFFFF.hpp"

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
compute_geom_1000(T&                               distributor,
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
        erirec::comp_electron_repulsion_geom1000_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_psdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_psff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ppsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ppdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ppff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_dsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_dpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    // leading [FS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_fsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fsdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fsff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ddsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ddpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ddff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_pfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_pfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_pfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_pfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_fpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fpsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fpdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fpff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_dfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_dfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_dfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_dfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_fdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_fdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_fddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_fdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [FF|XX] terms

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1000_ffss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_ffsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ffsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ffsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1000_ffpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ffpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ffpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1000_ffdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ffdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom1000_ffff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;

    std::abort();
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1000Func_hpp */
