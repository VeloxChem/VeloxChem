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

#ifndef NuclearPotentialGeom020SumRecGD_hpp
#define NuclearPotentialGeom020SumRecGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDD.hpp"
#include "NuclearPotentialGeom010PrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecFD.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialGeom020PrimRecDD.hpp"
#include "NuclearPotentialGeom020PrimRecDP.hpp"
#include "NuclearPotentialGeom020PrimRecDS.hpp"
#include "NuclearPotentialGeom020PrimRecFD.hpp"
#include "NuclearPotentialGeom020PrimRecFP.hpp"
#include "NuclearPotentialGeom020PrimRecGD.hpp"
#include "NuclearPotentialGeom020PrimRecPD.hpp"
#include "NuclearPotentialGeom020PrimRecPP.hpp"
#include "NuclearPotentialGeom020PrimRecPS.hpp"
#include "NuclearPotentialGeom020PrimRecSD.hpp"
#include "NuclearPotentialGeom020PrimRecSP.hpp"
#include "NuclearPotentialGeom020PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (G|AG(2)|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_gd(T&                               distributor,
                                       const CGtoBlock&                 bra_gto_block,
                                       const CGtoBlock&                 ket_gto_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices,
                                       const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external quadrupoles data

    const auto quadrupoles = distributor.data();

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

    CSimdArray<double> factors(20, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(4788, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(540, 1);

    CSimdArray<double> sbuffer(270, 1);

    // setup Boys function data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

        bf_data.set_active_width(ket_width);

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

                t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 9, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 9);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 9, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 12, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 15, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 18, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 21, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 24, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 27, 1, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 33, 2, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 39, 3, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 45, 4, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 51, 5, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 57, 6, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 63, 7, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 69, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 72, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 75, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 78, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 81, 2, 9, 12, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 90, 3, 12, 15, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 99, 4, 15, 18, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 108, 5, 18, 21, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 117, 6, 21, 24, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 126, 9, 27, 33, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 144, 12, 33, 39, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 162, 15, 39, 45, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 180, 18, 45, 51, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 198, 21, 51, 57, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 216, 24, 57, 63, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 234, 2, 3, 69, 72, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 240, 3, 4, 72, 75, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 246, 4, 5, 75, 78, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 252, 9, 12, 69, 81, 90, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 270, 12, 15, 72, 90, 99, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 288, 15, 18, 75, 99, 108, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 306, 18, 21, 78, 108, 117, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 324, 27, 33, 81, 126, 144, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 360, 33, 39, 90, 144, 162, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 396, 39, 45, 99, 162, 180, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 432, 45, 51, 108, 180, 198, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 468, 51, 57, 117, 198, 216, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 504, 2, 9, 12, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 513, 3, 12, 15, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 522, 4, 15, 18, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 531, 9, 27, 33, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 549, 12, 33, 39, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 567, 15, 39, 45, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 585, 18, 45, 51, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 603, 2, 3, 69, 72, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 612, 3, 4, 72, 75, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 621, 9, 12, 69, 81, 90, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 648, 12, 15, 72, 90, 99, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 675, 15, 18, 75, 99, 108, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 702, 27, 33, 81, 126, 144, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 756, 33, 39, 90, 144, 162, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 810, 39, 45, 99, 162, 180, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 864, 45, 51, 108, 180, 198, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 918, 69, 72, 234, 240, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 936, 72, 75, 240, 246, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 954, 81, 90, 234, 252, 270, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1008, 90, 99, 240, 270, 288, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1062, 99, 108, 246, 288, 306, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1116, 126, 144, 252, 324, 360, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1224, 144, 162, 270, 360, 396, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1332, 162, 180, 288, 396, 432, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1440, 180, 198, 306, 432, 468, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 1548, 27, 33, 504, 531, 549, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 1584, 33, 39, 513, 549, 567, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 1620, 39, 45, 522, 567, 585, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 1656, 81, 90, 504, 513, 603, 621, 648, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 1710, 90, 99, 513, 522, 612, 648, 675, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 1764, 126, 144, 531, 549, 621, 702, 756, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 1872, 144, 162, 549, 567, 648, 756, 810, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 1980, 162, 180, 567, 585, 675, 810, 864, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 2088, 234, 240, 603, 612, 918, 936, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 2124, 252, 270, 621, 648, 918, 954, 1008, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 2232, 270, 288, 648, 675, 936, 1008, 1062, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 2340, 324, 360, 702, 756, 954, 1116, 1224, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 2556, 360, 396, 756, 810, 1008, 1224, 1332, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 2772, 396, 432, 810, 864, 1062, 1332, 1440, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fp(pbuffer, 2988, 702, 756, 1548, 1584, 1656, 1764, 1872, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fp(pbuffer, 3168, 756, 810, 1584, 1620, 1710, 1872, 1980, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(pbuffer, 3348, 954, 1008, 1656, 1710, 2088, 2124, 2232, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fd(pbuffer, 3528, 1116, 1224, 1764, 1872, 2124, 2340, 2556, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fd(pbuffer, 3888, 1224, 1332, 1872, 1980, 2232, 2556, 2772, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_gd(pbuffer, 4248, 2340, 2556, 2988, 3168, 3348, 3528, 3888, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 4248, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom020SumRecGD_hpp */
