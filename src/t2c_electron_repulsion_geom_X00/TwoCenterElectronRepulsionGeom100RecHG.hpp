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

#ifndef TwoCenterElectronRepulsionGeom100RecHG_hpp
#define TwoCenterElectronRepulsionGeom100RecHG_hpp

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
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHF.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "TwoCenterElectronRepulsionPrimRecIG.hpp"
#include "GeometricalDerivatives1X0ForHY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)H|1/|r-r'||G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_hg(T& distributor,
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

    CSimdArray<double> pbuffer(5133, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(945, 1);

    CSimdArray<double> sbuffer(297, 1);

    // setup Boys function data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 0, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 9, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 10, bf_data, 10, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 11, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 14, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 17, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 20, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 23, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 26, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 29, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 32, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 35, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 38, 0, 1, 11, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 44, 1, 2, 14, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 50, 2, 3, 17, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 56, 3, 4, 20, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 62, 4, 5, 23, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 68, 5, 6, 26, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 74, 6, 7, 29, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 80, 7, 8, 32, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 86, 8, 9, 35, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 92, 11, 14, 50, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 102, 14, 17, 56, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 112, 17, 20, 62, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 122, 20, 23, 68, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 132, 23, 26, 74, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 142, 26, 29, 80, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 152, 29, 32, 86, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 162, 38, 44, 92, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 177, 44, 50, 102, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 192, 50, 56, 112, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 207, 56, 62, 122, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 222, 62, 68, 132, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 237, 68, 74, 142, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 252, 74, 80, 152, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 267, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 270, 4, 20, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 279, 5, 23, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 288, 6, 26, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 297, 20, 62, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 315, 23, 68, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 333, 26, 74, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 351, 50, 102, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 381, 56, 112, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 411, 62, 122, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 441, 68, 132, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 471, 74, 142, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 501, 102, 192, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 546, 112, 207, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 591, 122, 222, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 636, 132, 237, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 681, 142, 252, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 726, 4, 5, 267, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 732, 20, 23, 267, 288, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 750, 50, 56, 270, 297, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 786, 56, 62, 279, 315, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 822, 62, 68, 288, 333, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 858, 102, 112, 297, 411, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 918, 112, 122, 315, 441, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 978, 122, 132, 333, 471, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1038, 162, 177, 351, 501, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1128, 177, 192, 381, 546, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1218, 192, 207, 411, 591, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1308, 207, 222, 441, 636, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1398, 222, 237, 471, 681, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 1488, 270, 279, 726, 732, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1518, 297, 315, 732, 822, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1578, 351, 381, 750, 858, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1678, 381, 411, 786, 918, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1778, 411, 441, 822, 978, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 1878, 501, 546, 858, 1218, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2028, 546, 591, 918, 1308, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2178, 591, 636, 978, 1398, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 2328, 750, 786, 1488, 1518, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 2418, 858, 918, 1518, 1778, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 2568, 1038, 1128, 1578, 1878, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 2793, 1128, 1218, 1678, 2028, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 3018, 1218, 1308, 1778, 2178, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 3243, 1578, 1678, 2328, 2418, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 3453, 1878, 2028, 2418, 3018, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ig(pbuffer, 3768, 2568, 2793, 3243, 3453, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_hx(pbuffer, 4188, 2568, 3768, 1, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 4188, ket_width, ket_npgtos);
            }

            t2cfunc::transform<5, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 5, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecHG_hpp */
