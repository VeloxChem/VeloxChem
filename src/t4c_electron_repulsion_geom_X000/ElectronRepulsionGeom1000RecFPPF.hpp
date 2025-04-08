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

#ifndef ElectronRepulsionGeom1000RecFPPF_hpp
#define ElectronRepulsionGeom1000RecFPPF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
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

/// @brief Computes d^(1)/dA^(1)(FP|1/|r-r'||PF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1000_fppf(T& distributor,
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

    CSimdArray<double> pbuffer(3740, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1850, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2220, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(13524, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 10);

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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 9, pfactors, 16, bf_data, 9);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 155, 37, 43, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 170, 43, 49, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 185, 49, 55, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 200, 55, 61, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 61, 67, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 67, 73, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 245, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 248, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 251, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 260, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 269, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 278, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 296, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 314, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 332, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 350, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 380, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 410, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 440, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 470, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 500, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 545, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 590, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 635, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 680, 135, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 725, 3, 4, 245, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 731, 16, 19, 245, 251, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 749, 19, 22, 248, 260, 269, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 767, 43, 49, 251, 278, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 803, 49, 55, 260, 296, 314, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 839, 55, 61, 269, 314, 332, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 875, 85, 95, 278, 350, 380, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 935, 95, 105, 296, 380, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 995, 105, 115, 314, 410, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1055, 115, 125, 332, 440, 470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1115, 155, 170, 380, 500, 545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1205, 170, 185, 410, 545, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1295, 185, 200, 440, 590, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1385, 200, 215, 470, 635, 680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1475, 251, 260, 725, 731, 749, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1505, 278, 296, 731, 767, 803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1565, 296, 314, 749, 803, 839, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1625, 350, 380, 767, 875, 935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1725, 380, 410, 803, 935, 995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1825, 410, 440, 839, 995, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1925, 500, 545, 935, 1115, 1205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2075, 545, 590, 995, 1205, 1295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2225, 590, 635, 1055, 1295, 1385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2375, 767, 803, 1475, 1505, 1565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2465, 875, 935, 1505, 1625, 1725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2615, 935, 995, 1565, 1725, 1825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2765, 1115, 1205, 1725, 1925, 2075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2990, 1205, 1295, 1825, 2075, 2225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3215, 1625, 1725, 2375, 2465, 2615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 3425, 1925, 2075, 2615, 2765, 2990, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 350, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 500, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 875, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 1115, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 1625, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 325, pbuffer, 1925, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {350, 380});

                pbuffer.scale(2.0 * a_exp, {500, 545});

                pbuffer.scale(2.0 * a_exp, {875, 935});

                pbuffer.scale(2.0 * a_exp, {1115, 1205});

                pbuffer.scale(2.0 * a_exp, {1625, 1725});

                pbuffer.scale(2.0 * a_exp, {1925, 2075});

                pbuffer.scale(2.0 * a_exp, {2465, 2615});

                pbuffer.scale(2.0 * a_exp, {2765, 2990});

                pbuffer.scale(2.0 * a_exp, {3215, 3425});

                pbuffer.scale(2.0 * a_exp, {3425, 3740});

                t2cfunc::reduce(cbuffer, 475, pbuffer, 350, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 505, pbuffer, 500, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 875, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 1115, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 700, pbuffer, 1625, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 800, pbuffer, 1925, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 950, pbuffer, 2465, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1100, pbuffer, 2765, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1325, pbuffer, 3215, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1535, pbuffer, 3425, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 90, cbuffer, 75, 135, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 270, cbuffer, 225, 325, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 570, cbuffer, 475, 505, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 660, cbuffer, 550, 610, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 840, cbuffer, 700, 800, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1140, cbuffer, 950, 1100, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1590, cbuffer, 1325, 1535, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 0, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 63, ckbuffer, 90, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 189, ckbuffer, 270, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 1344, ckbuffer, 570, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 1407, ckbuffer, 660, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 1533, ckbuffer, 840, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 1743, ckbuffer, 1140, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 2058, ckbuffer, 1590, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 399, 0, 63, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 588, 63, 189, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 966, 399, 588, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2499, 1344, 1407, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2688, 1407, 1533, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3066, 1533, 1743, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3696, 1743, 2058, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 4641, 0, 2499, 2688, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 5208, 63, 2688, 3066, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 6342, 189, 3066, 3696, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 8232, 399, 4641, 5208, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 9366, 588, 5208, 6342, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_fpxx(skbuffer, 11634, 966, 8232, 9366, r_ab, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 11634, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 12264, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 882, skbuffer, 12894, 1, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 1, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1000RecFPPF_hpp */
