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

#ifndef ElectronRepulsionGeom1010RecDSPP_hpp
#define ElectronRepulsionGeom1010RecDSPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DS|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dspp(T& distributor,
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

    CSimdArray<double> pbuffer(696, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(528, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1296, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2997, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 95, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 98, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 101, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 110, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 119, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 128, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 146, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 164, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 182, 31, 55, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 212, 37, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 242, 43, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 272, 1, 2, 95, 98, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 278, 7, 10, 95, 101, 110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 296, 10, 13, 98, 110, 119, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 314, 25, 31, 110, 128, 146, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 350, 31, 37, 119, 146, 164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 386, 55, 65, 146, 182, 212, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 446, 65, 75, 164, 212, 242, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 506, 101, 110, 272, 278, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 536, 128, 146, 296, 314, 350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 596, 182, 212, 350, 386, 446, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 101, 9, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {101, 110});

                pbuffer.scale(2.0 * a_exp, {278, 296});

                pbuffer.scale(2.0 * a_exp, {506, 536});

                t2cfunc::reduce(cbuffer, 88, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 91, pbuffer, 101, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 278, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 118, pbuffer, 506, 30, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7, 10});

                pbuffer.scale(pfactors, 0, 2.0, {25, 31});

                pbuffer.scale(pfactors, 0, 2.0, {55, 65});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {101, 110});

                pbuffer.scale(pfactors, 0, 2.0, {128, 146});

                pbuffer.scale(pfactors, 0, 2.0, {182, 212});

                t2cfunc::reduce(cbuffer, 12, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 101, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 128, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 58, pbuffer, 182, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {55, 65});

                pbuffer.scale(2.0 * a_exp, {101, 110});

                pbuffer.scale(2.0 * a_exp, {128, 146});

                pbuffer.scale(2.0 * a_exp, {182, 212});

                pbuffer.scale(pfactors, 0, 2.0, {278, 296});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {314, 350});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {386, 446});

                pbuffer.scale(pfactors, 0, 2.0, {506, 536});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {536, 596});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {596, 696});

                t2cfunc::reduce(cbuffer, 148, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 151, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 157, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 167, pbuffer, 101, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 176, pbuffer, 128, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 182, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 224, pbuffer, 278, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 242, pbuffer, 314, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 278, pbuffer, 386, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 338, pbuffer, 506, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 368, pbuffer, 536, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 428, pbuffer, 596, 100, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 12, 15, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9, cbuffer, 15, 21, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 27, cbuffer, 0, 0, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 54, cbuffer, 31, 40, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 81, cbuffer, 40, 58, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 135, cbuffer, 3, 54, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 216, cbuffer, 148, 151, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 225, cbuffer, 151, 157, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 243, cbuffer, 88, 216, 225, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 270, cbuffer, 167, 176, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 297, cbuffer, 176, 194, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 351, cbuffer, 91, 270, 297, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 432, cbuffer, 224, 242, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 486, cbuffer, 242, 278, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 594, cbuffer, 100, 432, 486, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 756, cbuffer, 338, 368, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 846, cbuffer, 368, 428, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1026, cbuffer, 118, 756, 846, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 27, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 9, ckbuffer, 36, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 18, ckbuffer, 45, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 108, ckbuffer, 135, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 135, ckbuffer, 162, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 162, ckbuffer, 189, 0, 1);

            //t4cfunc::ket_transform<1, 1>(skbuffer, 918, ckbuffer, 0, 1, 0);

            //t4cfunc::ket_transform<1, 1>(skbuffer, 945, ckbuffer, 27, 1, 0);

            //t4cfunc::ket_transform<1, 1>(skbuffer, 972, ckbuffer, 54, 1, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2457, ckbuffer, 243, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2466, ckbuffer, 252, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2475, ckbuffer, 261, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2484, ckbuffer, 351, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2511, ckbuffer, 378, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2538, ckbuffer, 405, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2565, ckbuffer, 594, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2619, ckbuffer, 648, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2673, ckbuffer, 702, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2727, ckbuffer, 1026, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2817, ckbuffer, 1116, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2907, ckbuffer, 1206, 0, 3);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 918, 0, 108, r_ab, 1, 1);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 945, 9, 135, r_ab, 1, 1);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 972, 18, 162, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 27, 2457, 2484, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 189, 2484, 2565, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 432, 2565, 2727, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 999, 0, 27, 189, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1242, 108, 189, 432, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 1971, 918, 999, 1242, r_ab, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 1971, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 45, skbuffer, 2025, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 90, skbuffer, 2079, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 135, skbuffer, 2133, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 180, skbuffer, 2187, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 225, skbuffer, 2241, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 270, skbuffer, 2295, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 315, skbuffer, 2349, 1, 1);

            t4cfunc::bra_transform<2, 0>(sbuffer, 360, skbuffer, 2403, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDSPP_hpp */
