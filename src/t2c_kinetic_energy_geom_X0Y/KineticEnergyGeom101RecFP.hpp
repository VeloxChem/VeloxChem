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

#ifndef KineticEnergyGeom101RecFP_hpp
#define KineticEnergyGeom101RecFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "GeometricalDerivatives1X1ForFP.hpp"
#include "GtoBlock.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "KineticEnergyPrimRecDS.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "KineticEnergyPrimRecFP.hpp"
#include "KineticEnergyPrimRecFS.hpp"
#include "KineticEnergyPrimRecGD.hpp"
#include "KineticEnergyPrimRecGS.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecDS.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecFP.hpp"
#include "OverlapPrimRecFS.hpp"
#include "OverlapPrimRecGD.hpp"
#include "OverlapPrimRecGS.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)F|T|d^(1)/dB^(1)P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_geom_11_fp(T&                               distributor,
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

    CSimdArray<double> pbuffer(880, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(189, 1);

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

                ovlrec::comp_prim_overlap_sp(pbuffer, 2, 0, factors, 11);

                kinrec::comp_prim_kinetic_energy_sp(pbuffer, 5, 1, 2, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sd(pbuffer, 8, 0, 2, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sd(pbuffer, 14, 0, 1, 5, 8, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 20, 0, factors, 8);

                kinrec::comp_prim_kinetic_energy_ps(pbuffer, 23, 1, 20, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 26, 0, 2, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pp(pbuffer, 35, 1, 5, 26, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 44, 2, 8, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pd(pbuffer, 62, 5, 14, 44, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ds(pbuffer, 80, 0, 20, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_ds(pbuffer, 86, 0, 1, 23, 80, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 92, 2, 20, 26, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dp(pbuffer, 110, 2, 5, 23, 35, 92, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 128, 8, 26, 44, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dd(pbuffer, 164, 8, 14, 35, 62, 128, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fs(pbuffer, 200, 20, 80, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fs(pbuffer, 210, 20, 23, 86, 200, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fp(pbuffer, 220, 26, 80, 92, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fp(pbuffer, 250, 26, 35, 86, 110, 220, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 280, 44, 92, 128, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fd(pbuffer, 340, 44, 62, 110, 164, 280, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gs(pbuffer, 400, 80, 200, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gs(pbuffer, 415, 80, 86, 210, 400, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gd(pbuffer, 430, 128, 220, 280, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gd(pbuffer, 520, 128, 164, 250, 340, 430, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_11_fp(pbuffer, 610, 86, 164, 415, 520, 1, factors, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 610, ket_width, ket_npgtos);
            }

            t2cfunc::transform<3, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 1, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace kinrec

#endif /* KineticEnergyGeom101RecFP_hpp */
