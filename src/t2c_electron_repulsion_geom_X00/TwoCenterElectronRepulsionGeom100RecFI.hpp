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

#ifndef TwoCenterElectronRepulsionGeom100RecFI_hpp
#define TwoCenterElectronRepulsionGeom100RecFI_hpp

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
#include "TwoCenterElectronRepulsionPrimRecSI.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPI.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDI.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFI.hpp"
#include "TwoCenterElectronRepulsionPrimRecGI.hpp"
#include "GeometricalDerivatives1X0ForFY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)F|1/|r-r'||I)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_fi(T& distributor,
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

    CSimdArray<double> pbuffer(3498, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(840, 1);

    CSimdArray<double> sbuffer(273, 1);

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

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 267, 92, 102, 192, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 288, 102, 112, 207, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 309, 112, 122, 222, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 330, 122, 132, 237, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 351, 132, 142, 252, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 372, 162, 177, 267, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 400, 177, 192, 288, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 428, 192, 207, 309, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 456, 207, 222, 330, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 484, 222, 237, 351, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 512, 62, 122, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 542, 122, 222, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 587, 192, 288, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 650, 207, 309, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 713, 222, 330, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 776, 288, 428, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 860, 309, 456, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 944, 330, 484, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 1028, 192, 207, 512, 542, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 1118, 288, 309, 542, 713, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 1244, 372, 400, 587, 776, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 1412, 400, 428, 650, 860, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 1580, 428, 456, 713, 944, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 1748, 587, 650, 1028, 1118, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 1958, 776, 860, 1118, 1580, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 2238, 1244, 1412, 1748, 1958, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_fx(pbuffer, 2658, 1244, 2238, 1, 28, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 2658, ket_width, ket_npgtos);
            }

            t2cfunc::transform<3, 6>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 6, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecFI_hpp */
