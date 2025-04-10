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

#ifndef TwoCenterElectronRepulsionGeom100RecGD_hpp
#define TwoCenterElectronRepulsionGeom100RecGD_hpp

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
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"
#include "BoysFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)G|1/|r-r'||D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_gd(T& distributor,
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

    CSimdArray<double> pbuffer(1100, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(135, 1);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 7, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 7);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 7, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 10, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 25, 0, 1, 10, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 31, 1, 2, 13, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 37, 2, 3, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 43, 3, 4, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 49, 4, 5, 22, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 55, 2, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 58, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 61, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 64, 2, 13, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 73, 3, 16, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 82, 4, 19, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 91, 7, 25, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 109, 10, 31, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 127, 13, 37, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 145, 16, 43, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 163, 19, 49, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 181, 2, 3, 61, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 187, 7, 10, 55, 64, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 205, 10, 13, 58, 73, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 223, 13, 16, 61, 82, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 241, 25, 31, 64, 127, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 277, 31, 37, 73, 145, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 313, 37, 43, 82, 163, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 349, 55, 58, 181, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 359, 64, 73, 181, 223, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 389, 91, 109, 187, 241, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 449, 109, 127, 205, 277, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 509, 127, 145, 223, 313, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 569, 187, 205, 349, 359, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 614, 241, 277, 359, 509, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 704, 389, 449, 569, 614, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 830, 389, 704, 1, 6, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 830, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecGD_hpp */
