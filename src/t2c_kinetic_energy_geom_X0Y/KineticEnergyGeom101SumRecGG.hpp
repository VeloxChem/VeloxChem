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

#ifndef KineticEnergyGeom101SumRecGG_hpp
#define KineticEnergyGeom101SumRecGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "GeometricalDerivatives1X1ForGG.hpp"
#include "GtoBlock.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "KineticEnergyPrimRecDF.hpp"
#include "KineticEnergyPrimRecDG.hpp"
#include "KineticEnergyPrimRecDH.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "KineticEnergyPrimRecDS.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "KineticEnergyPrimRecFF.hpp"
#include "KineticEnergyPrimRecFG.hpp"
#include "KineticEnergyPrimRecFH.hpp"
#include "KineticEnergyPrimRecFP.hpp"
#include "KineticEnergyPrimRecGD.hpp"
#include "KineticEnergyPrimRecGF.hpp"
#include "KineticEnergyPrimRecGG.hpp"
#include "KineticEnergyPrimRecGH.hpp"
#include "KineticEnergyPrimRecHF.hpp"
#include "KineticEnergyPrimRecHH.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "KineticEnergyPrimRecPG.hpp"
#include "KineticEnergyPrimRecPH.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "KineticEnergyPrimRecSF.hpp"
#include "KineticEnergyPrimRecSG.hpp"
#include "KineticEnergyPrimRecSH.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "OverlapPrimRecDG.hpp"
#include "OverlapPrimRecDH.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecDS.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFG.hpp"
#include "OverlapPrimRecFH.hpp"
#include "OverlapPrimRecFP.hpp"
#include "OverlapPrimRecGD.hpp"
#include "OverlapPrimRecGF.hpp"
#include "OverlapPrimRecGG.hpp"
#include "OverlapPrimRecGH.hpp"
#include "OverlapPrimRecHF.hpp"
#include "OverlapPrimRecHH.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecPH.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "OverlapPrimRecSG.hpp"
#include "OverlapPrimRecSH.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)G|T|d^(1)/dB^(1)G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_kinetic_energy_geom_11_gg(T&                               distributor,
                                   const CGtoBlock&                 bra_gto_block,
                                   const CGtoBlock&                 ket_gto_block,
                                   const std::pair<size_t, size_t>& bra_indices,
                                   const std::pair<size_t, size_t>& ket_indices,
                                   const bool                       bra_eq_ket) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.coordinates();

    const auto ket_gto_exps = ket_gto_block.exponents();

    const auto ket_gto_norms = ket_gto_block.normalization_factors();

    const auto ket_gto_indices = ket_gto_block.orbital_indices();

    const auto ket_npgtos = ket_gto_block.number_of_primitives();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> factors(14, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(7107, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(2025, 1);

    CSimdArray<double> sbuffer(729, 1);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        factors.load(ket_gto_exps, ket_range, 0, ket_npgtos);

        factors.load(ket_gto_norms, ket_range, 1, ket_npgtos);

        factors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        sbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        pbuffer.set_active_width(ket_width);

        // loop over contracted basis functions on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            cbuffer.zero();

            sbuffer.zero();

            const auto r_a = bra_gto_coords[j];

            t2cfunc::comp_distances_ab(factors, 5, 2, r_a);

            for (size_t k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                kinrec::comp_prim_kinetic_energy_ss(pbuffer, 1, 0, factors, a_exp);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    ovlrec::comp_prim_overlap_sp(pbuffer, 2, 0, factors, 11);

                    kinrec::comp_prim_kinetic_energy_sp(pbuffer, 5, 1, 2, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sd(pbuffer, 8, 0, 2, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sd(pbuffer, 14, 0, 1, 5, 8, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sf(pbuffer, 20, 2, 8, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sf(pbuffer, 30, 2, 5, 14, 20, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sg(pbuffer, 40, 8, 20, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sg(pbuffer, 55, 8, 14, 30, 40, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sh(pbuffer, 70, 20, 40, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sh(pbuffer, 91, 20, 30, 55, 70, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_ps(pbuffer, 112, 0, factors, 8);

                    kinrec::comp_prim_kinetic_energy_ps(pbuffer, 115, 1, 112, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pp(pbuffer, 118, 0, 2, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pp(pbuffer, 127, 1, 5, 118, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pd(pbuffer, 136, 2, 8, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pd(pbuffer, 154, 5, 14, 136, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pf(pbuffer, 172, 8, 20, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pf(pbuffer, 202, 14, 30, 172, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pg(pbuffer, 232, 20, 40, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pg(pbuffer, 277, 30, 55, 232, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ph(pbuffer, 322, 40, 70, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ph(pbuffer, 385, 55, 91, 322, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ds(pbuffer, 448, 0, 112, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ds(pbuffer, 454, 0, 1, 115, 448, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dp(pbuffer, 460, 2, 112, 118, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dp(pbuffer, 478, 2, 5, 115, 127, 460, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dd(pbuffer, 496, 8, 118, 136, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dd(pbuffer, 532, 8, 14, 127, 154, 496, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_df(pbuffer, 568, 20, 136, 172, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_df(pbuffer, 628, 20, 30, 154, 202, 568, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dg(pbuffer, 688, 40, 172, 232, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dg(pbuffer, 778, 40, 55, 202, 277, 688, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dh(pbuffer, 868, 70, 232, 322, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dh(pbuffer, 994, 70, 91, 277, 385, 868, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fp(pbuffer, 1120, 118, 448, 460, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fp(pbuffer, 1150, 118, 127, 454, 478, 1120, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fd(pbuffer, 1180, 136, 460, 496, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fd(pbuffer, 1240, 136, 154, 478, 532, 1180, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ff(pbuffer, 1300, 172, 496, 568, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ff(pbuffer, 1400, 172, 202, 532, 628, 1300, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fg(pbuffer, 1500, 232, 568, 688, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fg(pbuffer, 1650, 232, 277, 628, 778, 1500, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fh(pbuffer, 1800, 322, 688, 868, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fh(pbuffer, 2010, 322, 385, 778, 994, 1800, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gd(pbuffer, 2220, 496, 1120, 1180, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gd(pbuffer, 2310, 496, 532, 1150, 1240, 2220, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gf(pbuffer, 2400, 568, 1180, 1300, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gf(pbuffer, 2550, 568, 628, 1240, 1400, 2400, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gg(pbuffer, 2700, 688, 1300, 1500, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gg(pbuffer, 2925, 688, 778, 1400, 1650, 2700, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gh(pbuffer, 3150, 868, 1500, 1800, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gh(pbuffer, 3465, 868, 994, 1650, 2010, 3150, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_hf(pbuffer, 3780, 1300, 2220, 2400, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_hf(pbuffer, 3990, 1300, 1400, 2310, 2550, 3780, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_hh(pbuffer, 4200, 1800, 2700, 3150, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_hh(pbuffer, 4641, 1800, 2010, 2925, 3465, 4200, factors, 8, a_exp);

                    t2cgeom::comp_prim_op_geom_11_gg(pbuffer, 5082, 1400, 2010, 3990, 4641, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace kinrec

#endif /* KineticEnergyGeom101SumRecGG_hpp */
