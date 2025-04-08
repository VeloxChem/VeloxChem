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

#ifndef ElectronRepulsionGeom1010RecPSDP_hpp
#define ElectronRepulsionGeom1010RecPSDP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(PS|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_psdp(T& distributor,
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

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(551, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(473, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2211, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1440, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(405, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            ckbuffer.zero();

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

                t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 7);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 95, 25, 31, 55, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 110, 31, 37, 65, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 125, 37, 43, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 140, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 143, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 152, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 161, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 179, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 197, 31, 55, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 227, 37, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 257, 65, 95, 110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 302, 75, 110, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 347, 7, 10, 140, 143, 152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 365, 25, 31, 152, 161, 179, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 401, 55, 65, 179, 197, 227, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 461, 95, 110, 227, 257, 302, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 25, 6, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {143, 152});

                pbuffer.scale(2.0 * a_exp, {161, 179});

                pbuffer.scale(2.0 * a_exp, {347, 365});

                pbuffer.scale(2.0 * a_exp, {365, 401});

                t2cfunc::reduce(cbuffer, 43, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 52, pbuffer, 143, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 61, pbuffer, 161, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 347, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 97, pbuffer, 365, 36, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7, 10});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {25, 31});

                pbuffer.scale(pfactors, 0, 2.0, {55, 65});

                pbuffer.scale(pfactors, 0, 2.0, {95, 110});

                t2cfunc::reduce(cbuffer, 9, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 28, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {55, 65});

                pbuffer.scale(2.0 * a_exp, {95, 110});

                pbuffer.scale(pfactors, 0, 2.0, {143, 152});

                pbuffer.scale(pfactors, 0, 2.0, {161, 179});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {197, 227});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {257, 302});

                pbuffer.scale(pfactors, 0, 2.0, {347, 365});

                pbuffer.scale(pfactors, 0, 2.0, {365, 401});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {401, 461});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {461, 551});

                t2cfunc::reduce(cbuffer, 133, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 136, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 142, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 152, pbuffer, 95, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 167, pbuffer, 143, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 176, pbuffer, 161, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 197, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 224, pbuffer, 257, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 269, pbuffer, 347, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 287, pbuffer, 365, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 323, pbuffer, 401, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 383, pbuffer, 461, 90, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 57, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 258, cbuffer, 43, 46, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 573, cbuffer, 52, 61, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1347, cbuffer, 79, 97, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 9, 12, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9, cbuffer, 12, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27, cbuffer, 18, 28, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 66, cbuffer, 0, 0, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 93, cbuffer, 3, 9, 27, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 147, 57, 66, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 201, cbuffer, 133, 136, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 210, cbuffer, 136, 142, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 228, cbuffer, 142, 152, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 267, cbuffer, 43, 201, 210, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 294, cbuffer, 46, 210, 228, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 348, 258, 267, 294, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 402, cbuffer, 167, 176, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 429, cbuffer, 176, 194, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 483, cbuffer, 194, 224, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 600, cbuffer, 52, 402, 429, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 681, cbuffer, 61, 429, 483, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 843, 573, 600, 681, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1005, cbuffer, 269, 287, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1059, cbuffer, 287, 323, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1167, cbuffer, 323, 383, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1401, cbuffer, 79, 1005, 1059, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1563, cbuffer, 97, 1059, 1167, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1887, 1347, 1401, 1563, cfactors, 6, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 147, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 15, ckbuffer, 165, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 30, ckbuffer, 183, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 990, ckbuffer, 348, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1005, ckbuffer, 366, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1020, ckbuffer, 384, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1035, ckbuffer, 843, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1080, ckbuffer, 897, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1125, ckbuffer, 951, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1170, ckbuffer, 1887, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1260, ckbuffer, 1995, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1350, ckbuffer, 2103, 0, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 45, 990, 1035, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 180, 1035, 1170, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 585, 0, 45, 180, r_ab, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 0, skbuffer, 585, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 45, skbuffer, 630, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 90, skbuffer, 675, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 135, skbuffer, 720, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 180, skbuffer, 765, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 225, skbuffer, 810, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 270, skbuffer, 855, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 315, skbuffer, 900, 2, 1);

            t4cfunc::bra_transform<1, 0>(sbuffer, 360, skbuffer, 945, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 0, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecPSDP_hpp */
