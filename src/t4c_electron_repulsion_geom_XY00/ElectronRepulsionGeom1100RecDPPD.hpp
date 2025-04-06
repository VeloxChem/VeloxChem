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

#ifndef ElectronRepulsionGeom1100RecDPPD_hpp
#define ElectronRepulsionGeom1100RecDPPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dppd(T& distributor,
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

    CSimdArray<double> pbuffer(2281, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1472, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3978, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(13170, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2025, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 9);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 8, pfactors, 16, bf_data, 8);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 135, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 138, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 141, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 144, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 153, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 162, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 171, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 180, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 198, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 216, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 234, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 252, 24, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 270, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 300, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 330, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 360, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 390, 63, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 420, 2, 3, 135, 138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 426, 3, 4, 138, 141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 432, 12, 15, 135, 144, 153, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 450, 15, 18, 138, 153, 162, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 468, 18, 21, 141, 162, 171, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 486, 33, 39, 144, 180, 198, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 522, 39, 45, 153, 198, 216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 558, 45, 51, 162, 216, 234, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 594, 51, 57, 171, 234, 252, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 630, 75, 85, 198, 270, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 690, 85, 95, 216, 300, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 750, 95, 105, 234, 330, 360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 810, 105, 115, 252, 360, 390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 870, 135, 138, 420, 426, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 880, 144, 153, 420, 432, 450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 910, 153, 162, 426, 450, 468, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 940, 180, 198, 432, 486, 522, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1000, 198, 216, 450, 522, 558, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1060, 216, 234, 468, 558, 594, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1120, 270, 300, 522, 630, 690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1220, 300, 330, 558, 690, 750, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1320, 330, 360, 594, 750, 810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1420, 432, 450, 870, 880, 910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1465, 486, 522, 880, 940, 1000, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1555, 522, 558, 910, 1000, 1060, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1645, 630, 690, 1000, 1120, 1220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1795, 690, 750, 1060, 1220, 1320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1945, 940, 1000, 1420, 1465, 1555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 2071, 1120, 1220, 1555, 1645, 1795, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 180, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 270, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {486, 522});

                pbuffer.scale(2.0 * b_exp, {630, 690});

                pbuffer.scale(2.0 * b_exp, {940, 1000});

                pbuffer.scale(2.0 * b_exp, {1120, 1220});

                t2cfunc::reduce(cbuffer, 64, pbuffer, 486, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 630, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 940, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1120, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {180, 198});

                pbuffer.scale(2.0 * a_exp, {270, 300});

                pbuffer.scale(a_exp / b_exp, {486, 522});

                pbuffer.scale(a_exp / b_exp, {630, 690});

                pbuffer.scale(a_exp / b_exp, {940, 1000});

                pbuffer.scale(a_exp / b_exp, {1120, 1220});

                t2cfunc::reduce(cbuffer, 320, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 326, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 336, pbuffer, 180, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 354, pbuffer, 270, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 384, pbuffer, 486, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 420, pbuffer, 630, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 480, pbuffer, 940, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 540, pbuffer, 1120, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {486, 522});

                pbuffer.scale(2.0 * b_exp, {630, 690});

                pbuffer.scale(2.0 * b_exp, {940, 1000});

                pbuffer.scale(2.0 * b_exp, {1120, 1220});

                pbuffer.scale(4.0 * a_exp * b_exp, {1465, 1555});

                pbuffer.scale(4.0 * a_exp * b_exp, {1645, 1795});

                pbuffer.scale(4.0 * a_exp * b_exp, {1945, 2071});

                pbuffer.scale(4.0 * a_exp * b_exp, {2071, 2281});

                t2cfunc::reduce(cbuffer, 640, pbuffer, 486, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 676, pbuffer, 630, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 940, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 796, pbuffer, 1120, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 896, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 986, pbuffer, 1645, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1136, pbuffer, 1945, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1262, pbuffer, 2071, 210, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18, cbuffer, 16, 34, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 558, cbuffer, 64, 100, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 666, cbuffer, 160, 220, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 846, cbuffer, 320, 326, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 864, cbuffer, 336, 354, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1080, cbuffer, 384, 420, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1512, cbuffer, 480, 540, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3042, cbuffer, 640, 676, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3150, cbuffer, 736, 796, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3330, cbuffer, 896, 986, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3600, cbuffer, 1136, 1262, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<1, 2>(skbuffer, 15, ckbuffer, 18, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 9915, ckbuffer, 558, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 10005, ckbuffer, 666, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 10155, ckbuffer, 846, 0, 0);

            t4cfunc::ket_transform<1, 2>(skbuffer, 10170, ckbuffer, 864, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 10350, ckbuffer, 1080, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 10710, ckbuffer, 1512, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 12390, ckbuffer, 3042, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 12480, ckbuffer, 3150, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 12630, ckbuffer, 3330, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 12855, ckbuffer, 3600, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 11985, 10170, 10350, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12120, 10350, 10710, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 3435, 15, 11985, 12120, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 60, 0, 9915, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 600, 15, 10005, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 3030, 15, 60, 600, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 10215, 10155, 12390, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 10440, 10170, 12480, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 10860, 10350, 12630, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 11310, 10710, 12855, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 195, 10170, 10215, 10440, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 870, 10350, 10440, 10860, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1680, 10710, 10860, 11310, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 3840, 60, 11985, 195, 870, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 5055, 600, 12120, 870, 1680, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 7485, 3030, 3435, 3840, 5055, r_ab, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 7485, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 225, skbuffer, 7755, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 450, skbuffer, 8025, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 675, skbuffer, 8295, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 900, skbuffer, 8565, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1125, skbuffer, 8835, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1350, skbuffer, 9105, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1575, skbuffer, 9375, 1, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1800, skbuffer, 9645, 1, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 1, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPPD_hpp */
