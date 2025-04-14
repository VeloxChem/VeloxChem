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

#ifndef TwoCenterElectronRepulsionRecID_hpp
#define TwoCenterElectronRepulsionRecID_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHP.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "TwoCenterElectronRepulsionPrimRecID.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (I|1/|r-r'||D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_id(T& distributor,
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

    CSimdArray<double> pbuffer(1435, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(168, 1);

    CSimdArray<double> sbuffer(65, 1);

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

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                t2cfunc::comp_boys_args_with_rho(bf_data, 9, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 9);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 0, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 30, 0, 1, 9, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 36, 1, 2, 12, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 42, 2, 3, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 48, 3, 4, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 54, 4, 5, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 60, 5, 6, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 66, 6, 7, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 72, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 75, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 78, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 81, 2, 12, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 90, 3, 15, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 99, 4, 18, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 108, 5, 21, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 117, 6, 24, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 126, 12, 42, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 144, 15, 48, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 162, 18, 54, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 180, 21, 60, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 198, 24, 66, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 216, 2, 3, 72, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 222, 3, 4, 75, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 228, 4, 5, 78, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 234, 12, 15, 72, 99, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 252, 15, 18, 75, 108, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 270, 18, 21, 78, 117, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 288, 30, 36, 81, 126, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 324, 36, 42, 90, 144, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 360, 42, 48, 99, 162, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 396, 48, 54, 108, 180, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 432, 54, 60, 117, 198, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 468, 72, 75, 228, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 478, 81, 90, 216, 234, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 508, 90, 99, 222, 252, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 538, 99, 108, 228, 270, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 568, 126, 144, 234, 360, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 628, 144, 162, 252, 396, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 688, 162, 180, 270, 432, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 748, 216, 222, 468, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 763, 234, 252, 468, 538, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 808, 288, 324, 478, 568, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 898, 324, 360, 508, 628, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 988, 360, 396, 538, 688, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 1078, 478, 508, 748, 763, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 1141, 568, 628, 763, 988, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_id(pbuffer, 1267, 808, 898, 1078, 1141, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 1267, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionRecID_hpp */
