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

#ifndef TwoCenterElectronRepulsionGeom100RecHH_hpp
#define TwoCenterElectronRepulsionGeom100RecHH_hpp

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
#include "TwoCenterElectronRepulsionPrimRecSH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHH.hpp"
#include "TwoCenterElectronRepulsionPrimRecIH.hpp"
#include "GeometricalDerivatives1X0ForHY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)H|1/|r-r'||H)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_hh(T& distributor,
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

    CSimdArray<double> pbuffer(7497, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1323, 1);

    CSimdArray<double> sbuffer(363, 1);

    // setup Boys function data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 11, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 11);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 9, bf_data, 10, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 10, bf_data, 11, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 11, 0, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 14, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 17, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 20, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 23, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 26, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 29, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 32, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 35, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 38, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 41, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 44, 0, 1, 17, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 50, 1, 2, 20, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 56, 2, 3, 23, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 62, 3, 4, 26, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 68, 4, 5, 29, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 74, 5, 6, 32, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 80, 6, 7, 35, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 86, 7, 8, 38, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 92, 8, 9, 41, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 98, 11, 14, 44, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 108, 14, 17, 50, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 118, 17, 20, 56, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 128, 20, 23, 62, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 138, 23, 26, 68, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 148, 26, 29, 74, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 158, 29, 32, 80, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 168, 32, 35, 86, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 178, 35, 38, 92, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 188, 44, 50, 118, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 203, 50, 56, 128, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 218, 56, 62, 138, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 233, 62, 68, 148, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 248, 68, 74, 158, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 263, 74, 80, 168, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 278, 80, 86, 178, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 293, 98, 108, 188, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 314, 108, 118, 203, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 335, 118, 128, 218, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 356, 128, 138, 233, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 377, 138, 148, 248, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 398, 148, 158, 263, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 419, 158, 168, 278, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 440, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 443, 5, 29, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 452, 23, 62, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 470, 26, 68, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 488, 29, 74, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 506, 62, 138, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 536, 68, 148, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 566, 74, 158, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 596, 118, 203, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 641, 128, 218, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 686, 138, 233, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 731, 148, 248, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 776, 158, 263, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 821, 203, 335, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 884, 218, 356, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 947, 233, 377, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1010, 248, 398, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1073, 263, 419, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 1136, 23, 26, 440, 443, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 1154, 62, 68, 443, 488, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1190, 118, 128, 452, 506, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1250, 128, 138, 470, 536, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1310, 138, 148, 488, 566, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1370, 203, 218, 506, 686, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1460, 218, 233, 536, 731, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1550, 233, 248, 566, 776, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1640, 293, 314, 596, 821, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1766, 314, 335, 641, 884, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1892, 335, 356, 686, 947, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2018, 356, 377, 731, 1010, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2144, 377, 398, 776, 1073, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 2270, 452, 470, 1136, 1154, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 2330, 506, 536, 1154, 1310, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2430, 596, 641, 1190, 1370, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2580, 641, 686, 1250, 1460, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 2730, 686, 731, 1310, 1550, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 2880, 821, 884, 1370, 1892, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3090, 884, 947, 1460, 2018, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 3300, 947, 1010, 1550, 2144, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 3510, 1190, 1250, 2270, 2330, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 3660, 1370, 1460, 2330, 2730, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 3885, 1640, 1766, 2430, 2880, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 4200, 1766, 1892, 2580, 3090, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 4515, 1892, 2018, 2730, 3300, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 4830, 2430, 2580, 3510, 3660, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 5145, 2880, 3090, 3660, 4515, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ih(pbuffer, 5586, 3885, 4200, 4830, 5145, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_hx(pbuffer, 6174, 3885, 5586, 1, 21, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 6174, ket_width, ket_npgtos);
            }

            t2cfunc::transform<5, 5>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 5, 5, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecHH_hpp */
