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

#ifndef TwoCenterElectronRepulsionRecIP_hpp
#define TwoCenterElectronRepulsionRecIP_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecHS.hpp"
#include "TwoCenterElectronRepulsionPrimRecHP.hpp"
#include "TwoCenterElectronRepulsionPrimRecIP.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (I|1/|r-r'||P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_ip(T& distributor,
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

    CSimdArray<double> pbuffer(634, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(84, 1);

    CSimdArray<double> sbuffer(39, 1);

    // setup Boys function data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 8, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 8);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 7, 0, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 10, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 25, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 28, 1, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 31, 2, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 34, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 37, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 40, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 43, 1, 13, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 52, 2, 16, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 61, 3, 19, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 70, 4, 22, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 79, 5, 25, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 88, 1, 2, 34, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 94, 2, 3, 37, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 100, 3, 4, 40, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 106, 7, 10, 28, 43, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 124, 10, 13, 31, 52, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 142, 13, 16, 34, 61, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 160, 16, 19, 37, 70, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 178, 19, 22, 40, 79, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 196, 28, 31, 88, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 206, 31, 34, 94, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 216, 34, 37, 100, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 226, 43, 52, 88, 142, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 256, 52, 61, 94, 160, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 286, 61, 70, 100, 178, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 316, 88, 94, 216, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 331, 106, 124, 196, 226, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 376, 124, 142, 206, 256, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 421, 142, 160, 216, 286, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hs(pbuffer, 466, 196, 206, 316, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 487, 226, 256, 316, 421, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ip(pbuffer, 550, 331, 376, 466, 487, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 550, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 1, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionRecIP_hpp */
