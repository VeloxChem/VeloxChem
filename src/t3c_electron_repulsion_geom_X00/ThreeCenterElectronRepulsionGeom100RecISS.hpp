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

#ifndef ThreeCenterElectronRepulsionGeom100RecISS_hpp
#define ThreeCenterElectronRepulsionGeom100RecISS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecIXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T3CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"
#include "GtoBlock.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dA^(1)(I|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_iss(T& distributor,
                                    const CGtoBlock& bra_gto_block,
                                    const CGtoPairBlock& ket_gto_pair_block,
                                    const std::pair<size_t, size_t>& bra_range) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto c_coords = ket_gto_pair_block.bra_coordinates();

    const auto d_coords = ket_gto_pair_block.ket_coordinates();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(23, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(280, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(177, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(39, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(39, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_gto_pair_block.number_of_contracted_pairs();

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), size_t{0});

        pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);

        pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);

        pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);

        pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);

        pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);

        pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_range.first; j < bra_range.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = bra_gto_coords[j];

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

                t3cfunc::comp_distances_aq(pfactors, 13, 10, r_a);

                t3cfunc::comp_coordinates_w(pfactors, 17, 10, r_a, a_exp);

                t4cfunc::comp_distances_wp(pfactors, 20, 17, r_a);

                t3cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 8);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 7, 0, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 10, 1, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 13, 2, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 16, 3, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 19, 4, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 22, 5, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 25, 6, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 28, 0, 1, 13, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 34, 1, 2, 16, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 40, 2, 3, 19, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 46, 3, 4, 22, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 52, 4, 5, 25, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 58, 7, 10, 28, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 68, 10, 13, 34, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 78, 13, 16, 40, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 88, 16, 19, 46, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 98, 19, 22, 52, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 108, 28, 34, 78, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 123, 34, 40, 88, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 138, 40, 46, 98, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 153, 58, 68, 108, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 174, 68, 78, 123, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 195, 78, 88, 138, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_iss(pbuffer, 216, 108, 123, 195, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_kss(pbuffer, 244, 153, 174, 216, pfactors, 20, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 153, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 105, pbuffer, 244, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {244, 280});

                t2cfunc::reduce(cbuffer, 141, pbuffer, 244, 36, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 21, 0, 141, 0, 0);

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 21, 0, 0);

            t3cfunc::bra_transform<6>(skbuffer, 13, cbuffer, 49, 0, 0);

            t3cfunc::bra_transform<6>(skbuffer, 26, cbuffer, 77, 0, 0);

            t3cfunc::ket_transform<0, 0>(sbuffer, 0, skbuffer, 0, 6);

            t3cfunc::ket_transform<0, 0>(sbuffer, 13, skbuffer, 13, 6);

            t3cfunc::ket_transform<0, 0>(sbuffer, 26, skbuffer, 26, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 0, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecISS_hpp */
