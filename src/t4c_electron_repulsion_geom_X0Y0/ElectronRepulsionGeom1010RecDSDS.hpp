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

#ifndef ElectronRepulsionGeom1010RecDSDS_hpp
#define ElectronRepulsionGeom1010RecDSDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
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
#include "ElectronRepulsionPrimRecSFSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DS|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dsds(T& distributor,
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

    CSimdArray<double> pbuffer(715, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(576, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2088, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1665, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(225, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 95, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 98, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 101, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 104, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 113, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 122, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 131, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 149, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 167, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 185, 31, 55, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 215, 37, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 245, 43, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 275, 0, 1, 95, 98, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 281, 1, 2, 98, 101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 287, 7, 10, 98, 104, 113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 305, 10, 13, 101, 113, 122, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 323, 25, 31, 113, 131, 149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 359, 31, 37, 122, 149, 167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 395, 55, 65, 149, 185, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 455, 65, 75, 167, 215, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 515, 95, 98, 275, 281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 525, 104, 113, 281, 287, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 555, 131, 149, 305, 323, 359, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 615, 185, 215, 359, 395, 455, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 95, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7, pbuffer, 104, 9, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {95, 98});

                pbuffer.scale(2.0 * a_exp, {104, 113});

                pbuffer.scale(2.0 * a_exp, {275, 281});

                pbuffer.scale(2.0 * a_exp, {287, 305});

                pbuffer.scale(2.0 * a_exp, {515, 525});

                pbuffer.scale(2.0 * a_exp, {525, 555});

                t2cfunc::reduce(cbuffer, 96, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 97, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 95, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 103, pbuffer, 104, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 112, pbuffer, 275, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 118, pbuffer, 287, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 136, pbuffer, 515, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 146, pbuffer, 525, 30, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {0, 1});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7, 10});

                pbuffer.scale(pfactors, 0, 2.0, {25, 31});

                pbuffer.scale(pfactors, 0, 2.0, {55, 65});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {95, 98});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {104, 113});

                pbuffer.scale(pfactors, 0, 2.0, {131, 149});

                pbuffer.scale(pfactors, 0, 2.0, {185, 215});

                t2cfunc::reduce(cbuffer, 16, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 26, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 95, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 39, pbuffer, 104, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 131, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 66, pbuffer, 185, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {55, 65});

                pbuffer.scale(2.0 * a_exp, {95, 98});

                pbuffer.scale(2.0 * a_exp, {104, 113});

                pbuffer.scale(2.0 * a_exp, {131, 149});

                pbuffer.scale(2.0 * a_exp, {185, 215});

                pbuffer.scale(pfactors, 0, 2.0, {275, 281});

                pbuffer.scale(pfactors, 0, 2.0, {287, 305});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {323, 359});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {395, 455});

                pbuffer.scale(pfactors, 0, 2.0, {515, 525});

                pbuffer.scale(pfactors, 0, 2.0, {525, 555});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {555, 615});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {615, 715});

                t2cfunc::reduce(cbuffer, 176, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 177, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 180, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 196, pbuffer, 95, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 199, pbuffer, 104, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 208, pbuffer, 131, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 226, pbuffer, 185, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 256, pbuffer, 275, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 262, pbuffer, 287, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 323, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 316, pbuffer, 395, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 376, pbuffer, 515, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 386, pbuffer, 525, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 416, pbuffer, 555, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 476, pbuffer, 615, 100, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 30, cbuffer, 0, 1, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 177, cbuffer, 4, 7, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 378, cbuffer, 96, 97, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 525, cbuffer, 100, 103, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 876, cbuffer, 112, 118, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1518, cbuffer, 136, 146, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 16, 17, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3, cbuffer, 17, 20, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12, cbuffer, 20, 26, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 33, cbuffer, 0, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 42, cbuffer, 1, 3, 12, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 69, 30, 33, 42, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 87, cbuffer, 36, 39, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 96, cbuffer, 39, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 123, cbuffer, 48, 66, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 186, cbuffer, 4, 87, 96, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 213, cbuffer, 7, 96, 123, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 294, 177, 186, 213, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 348, cbuffer, 176, 177, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 351, cbuffer, 177, 180, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 360, cbuffer, 180, 186, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 381, cbuffer, 96, 348, 351, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 390, cbuffer, 97, 351, 360, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 417, 378, 381, 390, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 435, cbuffer, 196, 199, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 444, cbuffer, 199, 208, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 471, cbuffer, 208, 226, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 534, cbuffer, 100, 435, 444, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 561, cbuffer, 103, 444, 471, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 642, 525, 534, 561, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 696, cbuffer, 256, 262, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 714, cbuffer, 262, 280, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 768, cbuffer, 280, 316, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 894, cbuffer, 112, 696, 714, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 948, cbuffer, 118, 714, 768, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1110, 876, 894, 948, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1218, cbuffer, 376, 386, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1248, cbuffer, 386, 416, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1338, cbuffer, 416, 476, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1548, cbuffer, 136, 1218, 1248, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1638, cbuffer, 146, 1248, 1338, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1908, 1518, 1548, 1638, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 69, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 5, ckbuffer, 75, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 10, ckbuffer, 81, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 60, ckbuffer, 294, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 75, ckbuffer, 312, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 90, ckbuffer, 330, 0, 1);

            //t4cfunc::ket_transform<2, 0>(skbuffer, 510, ckbuffer, 0, 1, 0);

            //t4cfunc::ket_transform<2, 0>(skbuffer, 525, ckbuffer, 18, 1, 0);

            //t4cfunc::ket_transform<2, 0>(skbuffer, 540, ckbuffer, 36, 1, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1365, ckbuffer, 417, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1370, ckbuffer, 423, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1375, ckbuffer, 429, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1380, ckbuffer, 642, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1395, ckbuffer, 660, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1410, ckbuffer, 678, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1425, ckbuffer, 1110, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1455, ckbuffer, 1146, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1485, ckbuffer, 1182, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1515, ckbuffer, 1908, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1565, ckbuffer, 1968, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1615, ckbuffer, 2028, 0, 3);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 510, 0, 60, r_ab, 2, 0);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 525, 5, 75, r_ab, 2, 0);
            
            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 540, 10, 90, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 15, 1365, 1380, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 105, 1380, 1425, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 240, 1425, 1515, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 555, 0, 15, 105, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 690, 60, 105, 240, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 1095, 510, 555, 690, r_ab, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 1095, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 25, skbuffer, 1125, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 50, skbuffer, 1155, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 75, skbuffer, 1185, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 100, skbuffer, 1215, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 125, skbuffer, 1245, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 150, skbuffer, 1275, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 175, skbuffer, 1305, 2, 0);

            t4cfunc::bra_transform<2, 0>(sbuffer, 200, skbuffer, 1335, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDSDS_hpp */
