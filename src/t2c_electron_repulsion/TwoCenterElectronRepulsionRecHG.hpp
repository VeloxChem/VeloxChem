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

#ifndef TwoCenterElectronRepulsionRecHG_hpp
#define TwoCenterElectronRepulsionRecHG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecSF.hpp"
#include "TwoCenterElectronRepulsionPrimRecSG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (H|1/|r-r'||G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_hg(T& distributor,
                           const CGtoBlock& bra_gto_block,
                           const CGtoBlock& ket_gto_block,
                           const std::pair<size_t, size_t>& bra_indices,
                           const std::pair<size_t, size_t>& ket_indices,
                           const bool bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(2395, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(315, 1);

    CSimdArray<double> sbuffer(99, 1);

    // setup Boys function data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                t2cfunc::comp_boys_args_with_rho(bf_data, 10, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 10);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 30, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 33, 0, 1, 12, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 39, 1, 2, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 45, 2, 3, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 51, 3, 4, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 57, 4, 5, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 63, 5, 6, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 69, 6, 7, 30, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 75, 9, 12, 39, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 85, 12, 15, 45, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 95, 15, 18, 51, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 105, 18, 21, 57, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 115, 21, 24, 63, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 125, 24, 27, 69, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 135, 33, 39, 85, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 150, 39, 45, 95, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 165, 45, 51, 105, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 180, 51, 57, 115, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 195, 57, 63, 125, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 210, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 213, 4, 21, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 222, 15, 45, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 240, 18, 51, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 258, 21, 57, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 276, 45, 95, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 306, 51, 105, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 336, 57, 115, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 366, 75, 135, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 411, 85, 150, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 456, 95, 165, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 501, 105, 180, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 546, 115, 195, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 591, 15, 18, 210, 213, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 609, 45, 51, 213, 258, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 645, 75, 85, 222, 276, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 705, 85, 95, 240, 306, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 765, 95, 105, 258, 336, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 825, 135, 150, 276, 456, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 915, 150, 165, 306, 501, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1005, 165, 180, 336, 546, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1095, 222, 240, 591, 609, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1155, 276, 306, 609, 765, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 1255, 366, 411, 645, 825, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 1405, 411, 456, 705, 915, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 1555, 456, 501, 765, 1005, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 1705, 645, 705, 1095, 1155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 1855, 825, 915, 1155, 1555, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 2080, 1255, 1405, 1705, 1855, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 2080, ket_width, ket_npgtos);
            }

            t2cfunc::transform<5, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 5, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionRecHG_hpp */
