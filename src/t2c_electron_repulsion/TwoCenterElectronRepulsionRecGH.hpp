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

#ifndef TwoCenterElectronRepulsionRecGH_hpp
#define TwoCenterElectronRepulsionRecGH_hpp

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
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (G|1/|r-r'||H)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_gh(T& distributor,
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

    CSimdArray<double> pbuffer(1903, ket_npgtos);

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

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 0, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 30, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 33, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 36, 0, 1, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 42, 1, 2, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 48, 2, 3, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 54, 3, 4, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 60, 4, 5, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 66, 5, 6, 30, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 72, 6, 7, 33, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 78, 9, 12, 36, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 88, 12, 15, 42, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 98, 15, 18, 48, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 108, 18, 21, 54, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 118, 21, 24, 60, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 128, 24, 27, 66, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 138, 27, 30, 72, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 148, 36, 42, 98, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 163, 42, 48, 108, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 178, 48, 54, 118, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 193, 54, 60, 128, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 208, 60, 66, 138, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 223, 78, 88, 148, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 244, 88, 98, 163, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 265, 98, 108, 178, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 286, 108, 118, 193, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 307, 118, 128, 208, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 328, 21, 54, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 346, 54, 118, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 376, 98, 163, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 421, 108, 178, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 466, 118, 193, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 511, 163, 265, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 574, 178, 286, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 637, 193, 307, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 700, 98, 108, 328, 346, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 760, 163, 178, 346, 466, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 850, 223, 244, 376, 511, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 976, 244, 265, 421, 574, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1102, 265, 286, 466, 637, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 1228, 376, 421, 700, 760, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 1378, 511, 574, 760, 1102, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 1588, 850, 976, 1228, 1378, factors, 8, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 1588, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 5>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 5, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionRecGH_hpp */
