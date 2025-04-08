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

#ifndef KineticEnergyGeom101SumRecDG_hpp
#define KineticEnergyGeom101SumRecDG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "GeometricalDerivatives1X1ForDG.hpp"
#include "GtoBlock.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "KineticEnergyPrimRecDF.hpp"
#include "KineticEnergyPrimRecDG.hpp"
#include "KineticEnergyPrimRecDH.hpp"
#include "KineticEnergyPrimRecFF.hpp"
#include "KineticEnergyPrimRecFH.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "KineticEnergyPrimRecPG.hpp"
#include "KineticEnergyPrimRecPH.hpp"
#include "KineticEnergyPrimRecPP.hpp"
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
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFH.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecPH.hpp"
#include "OverlapPrimRecPP.hpp"
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

/// @brief Computes (d^(1)/dA^(1)D|T|d^(1)/dB^(1)G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_kinetic_energy_geom_11_dg(T&                               distributor,
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

    CSimdArray<double> pbuffer(2496, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(810, 1);

    CSimdArray<double> sbuffer(405, 1);

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

                    ovlrec::comp_prim_overlap_pp(pbuffer, 112, 0, 2, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pp(pbuffer, 121, 1, 5, 112, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pd(pbuffer, 130, 2, 8, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pd(pbuffer, 148, 5, 14, 130, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pf(pbuffer, 166, 8, 20, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pf(pbuffer, 196, 14, 30, 166, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pg(pbuffer, 226, 20, 40, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pg(pbuffer, 271, 30, 55, 226, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ph(pbuffer, 316, 40, 70, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ph(pbuffer, 379, 55, 91, 316, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dd(pbuffer, 442, 8, 112, 130, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dd(pbuffer, 478, 8, 14, 121, 148, 442, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_df(pbuffer, 514, 20, 130, 166, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_df(pbuffer, 574, 20, 30, 148, 196, 514, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dg(pbuffer, 634, 40, 166, 226, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dg(pbuffer, 724, 40, 55, 196, 271, 634, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dh(pbuffer, 814, 70, 226, 316, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dh(pbuffer, 940, 70, 91, 271, 379, 814, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ff(pbuffer, 1066, 166, 442, 514, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ff(pbuffer, 1166, 166, 196, 478, 574, 1066, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fh(pbuffer, 1266, 316, 634, 814, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fh(pbuffer, 1476, 316, 379, 724, 940, 1266, factors, 8, a_exp);

                    t2cgeom::comp_prim_op_geom_11_dg(pbuffer, 1686, 196, 379, 1166, 1476, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<2, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace kinrec

#endif /* KineticEnergyGeom101SumRecDG_hpp */
