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

#ifndef ElectronRepulsionGeom1100RecDSSP_hpp
#define ElectronRepulsionGeom1100RecDSSP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DS|1/|r-r'||SP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dssp(T& distributor,
                                      const CGtoPairBlock& bra_gto_pair_block,
                                      const CGtoPairBlock& ket_gto_pair_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords = bra_gto_pair_block.bra_coordinates();

    const auto b_coords = bra_gto_pair_block.ket_coordinates();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

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

    CSimdArray<double> pfactors(29, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(247, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(162, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1206, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(135, 1);

    // setup Boys fuction data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

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

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = a_coords[j];

            const auto r_b = b_coords[j];

            const auto a_xyz = r_a.coordinates();

            const auto b_xyz = r_b.coordinates();

            const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = a_vec_exps[k * bra_ncgtos + j];

                const auto b_exp = b_vec_exps[k * bra_ncgtos + j];

                const auto ab_norm = ab_vec_norms[k * bra_ncgtos + j];

                const auto ab_ovl = ab_vec_ovls[k * bra_ncgtos + j];

                const auto p_x = (a_xyz[0] * a_exp + b_xyz[0] * b_exp) / (a_exp + b_exp);

                const auto p_y = (a_xyz[1] * a_exp + b_xyz[1] * b_exp) / (a_exp + b_exp);

                const auto p_z = (a_xyz[2] * a_exp + b_xyz[2] * b_exp) / (a_exp + b_exp);

                const auto r_p = TPoint<double>({p_x, p_y, p_z});

                const auto pb_x = p_x - b_xyz[0];

                const auto pb_y = p_y - b_xyz[1];

                const auto pb_z = p_z - b_xyz[2];

                const auto r_pb = TPoint<double>({pb_x, pb_y, pb_z});

                t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

                t4cfunc::comp_distances_pq(pfactors, 13, 10, r_p);

                t4cfunc::comp_coordinates_w(pfactors, 17, 10, r_p, a_exp, b_exp);

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

                t4cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 6);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 6, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 21, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 24, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 27, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 30, 1, 6, 9, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 39, 2, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 48, 3, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 57, 4, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 66, 1, 2, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 72, 2, 3, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 78, 6, 9, 21, 30, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 96, 9, 12, 24, 39, 48, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 114, 12, 15, 27, 48, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 132, 21, 24, 66, 72, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 142, 30, 39, 66, 78, 96, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 172, 39, 48, 72, 96, 114, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 202, 78, 96, 132, 142, 172, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 6, 3, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {30, 39});

                pbuffer.scale(2.0 * b_exp, {78, 96});

                t2cfunc::reduce(cbuffer, 3, pbuffer, 30, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 78, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(a_exp / b_exp, {30, 39});

                pbuffer.scale(a_exp / b_exp, {78, 96});

                t2cfunc::reduce(cbuffer, 30, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 33, pbuffer, 30, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 78, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {30, 39});

                pbuffer.scale(2.0 * b_exp, {78, 96});

                pbuffer.scale(4.0 * a_exp * b_exp, {142, 172});

                pbuffer.scale(4.0 * a_exp * b_exp, {202, 247});

                t2cfunc::reduce(cbuffer, 60, pbuffer, 30, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 69, pbuffer, 78, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 87, pbuffer, 142, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 117, pbuffer, 202, 45, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 1>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 1>(skbuffer, 813, cbuffer, 3, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 822, cbuffer, 12, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 840, cbuffer, 30, 0, 0);

            t4cfunc::ket_transform<0, 1>(skbuffer, 843, cbuffer, 33, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 879, cbuffer, 42, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 1077, cbuffer, 60, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 1086, cbuffer, 69, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 1104, cbuffer, 87, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 1134, cbuffer, 117, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1041, 840, 843, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1050, 843, 879, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1179, 1077, 1086, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 300, 0, 1041, 1050, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 3, 0, 822, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_psxx(skbuffer, 273, 0, 813, 3, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 852, 840, 1086, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 897, 843, 1104, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 951, 879, 1134, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 30, 843, 852, 897, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 111, 879, 897, 951, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_psxx(skbuffer, 327, 813, 1041, 1179, 30, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 408, 3, 1050, 30, 111, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dsxx(skbuffer, 651, 273, 300, 327, 408, r_ab, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 651, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 15, skbuffer, 669, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 30, skbuffer, 687, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 45, skbuffer, 705, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 60, skbuffer, 723, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 75, skbuffer, 741, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 90, skbuffer, 759, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 105, skbuffer, 777, 0, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 120, skbuffer, 795, 0, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 0, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDSSP_hpp */
