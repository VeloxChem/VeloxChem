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

#ifndef ElectronRepulsionGeom1000RecFDDD_hpp
#define ElectronRepulsionGeom1000RecFDDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
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

/// @brief Computes d^(1)/dA^(1)(FD|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1000_fddd(T& distributor,
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

    CSimdArray<double> pbuffer(6998, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3441, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(9324, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(27450, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2625, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 11);

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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 10, pfactors, 16, bf_data, 10);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 280, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 283, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 286, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 289, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 292, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 301, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 310, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 319, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 328, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 337, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 355, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 373, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 391, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 409, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 427, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 445, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 475, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 505, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 535, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 565, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 595, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 625, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 670, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 715, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 760, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 805, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 850, 155, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 895, 2, 3, 280, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 901, 3, 4, 283, 286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 907, 4, 5, 286, 289, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 913, 14, 17, 280, 292, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 931, 17, 20, 283, 301, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 949, 20, 23, 286, 310, 319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 967, 23, 26, 289, 319, 328, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 985, 41, 47, 292, 337, 355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1021, 47, 53, 301, 355, 373, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1057, 53, 59, 310, 373, 391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1093, 59, 65, 319, 391, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1129, 65, 71, 328, 409, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1165, 95, 105, 355, 445, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1225, 105, 115, 373, 475, 505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1285, 115, 125, 391, 505, 535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1345, 125, 135, 409, 535, 565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1405, 135, 145, 427, 565, 595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1465, 175, 190, 475, 625, 670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1555, 190, 205, 505, 670, 715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1645, 205, 220, 535, 715, 760, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1735, 220, 235, 565, 760, 805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1825, 235, 250, 595, 805, 850, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1915, 280, 283, 895, 901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1925, 283, 286, 901, 907, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1935, 292, 301, 895, 913, 931, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1965, 301, 310, 901, 931, 949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1995, 310, 319, 907, 949, 967, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2025, 337, 355, 913, 985, 1021, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2085, 355, 373, 931, 1021, 1057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2145, 373, 391, 949, 1057, 1093, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2205, 391, 409, 967, 1093, 1129, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2265, 445, 475, 1021, 1165, 1225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2365, 475, 505, 1057, 1225, 1285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2465, 505, 535, 1093, 1285, 1345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2565, 535, 565, 1129, 1345, 1405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2665, 625, 670, 1225, 1465, 1555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2815, 670, 715, 1285, 1555, 1645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2965, 715, 760, 1345, 1645, 1735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3115, 760, 805, 1405, 1735, 1825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3265, 895, 901, 1915, 1925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3280, 913, 931, 1915, 1935, 1965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3325, 931, 949, 1925, 1965, 1995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3370, 985, 1021, 1935, 2025, 2085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3460, 1021, 1057, 1965, 2085, 2145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3550, 1057, 1093, 1995, 2145, 2205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3640, 1165, 1225, 2085, 2265, 2365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3790, 1225, 1285, 2145, 2365, 2465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3940, 1285, 1345, 2205, 2465, 2565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4090, 1465, 1555, 2365, 2665, 2815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4315, 1555, 1645, 2465, 2815, 2965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4540, 1645, 1735, 2565, 2965, 3115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4765, 1935, 1965, 3265, 3280, 3325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4828, 2025, 2085, 3280, 3370, 3460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4954, 2085, 2145, 3325, 3460, 3550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5080, 2265, 2365, 3460, 3640, 3790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5290, 2365, 2465, 3550, 3790, 3940, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5500, 2665, 2815, 3790, 4090, 4315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5815, 2815, 2965, 3940, 4315, 4540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 6130, 3370, 3460, 4765, 4828, 4954, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6298, 3640, 3790, 4954, 5080, 5290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 6578, 4090, 4315, 5290, 5500, 5815, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 2025, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 246, pbuffer, 2265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 346, pbuffer, 2665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 3370, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 586, pbuffer, 3640, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 4090, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {985, 1021});

                pbuffer.scale(2.0 * a_exp, {1165, 1225});

                pbuffer.scale(2.0 * a_exp, {1465, 1555});

                pbuffer.scale(2.0 * a_exp, {2025, 2085});

                pbuffer.scale(2.0 * a_exp, {2265, 2365});

                pbuffer.scale(2.0 * a_exp, {2665, 2815});

                pbuffer.scale(2.0 * a_exp, {3370, 3460});

                pbuffer.scale(2.0 * a_exp, {3640, 3790});

                pbuffer.scale(2.0 * a_exp, {4090, 4315});

                pbuffer.scale(2.0 * a_exp, {4828, 4954});

                pbuffer.scale(2.0 * a_exp, {5080, 5290});

                pbuffer.scale(2.0 * a_exp, {5500, 5815});

                pbuffer.scale(2.0 * a_exp, {6130, 6298});

                pbuffer.scale(2.0 * a_exp, {6298, 6578});

                pbuffer.scale(2.0 * a_exp, {6578, 6998});

                t2cfunc::reduce(cbuffer, 961, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 997, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1057, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1147, pbuffer, 2025, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1207, pbuffer, 2265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1307, pbuffer, 2665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1457, pbuffer, 3370, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1547, pbuffer, 3640, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1697, pbuffer, 4090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1922, pbuffer, 4828, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2048, pbuffer, 5080, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2258, pbuffer, 5500, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2573, pbuffer, 6130, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2741, pbuffer, 6298, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3021, pbuffer, 6578, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 108, cbuffer, 36, 96, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 288, 0, 108, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 504, cbuffer, 186, 246, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 684, cbuffer, 246, 346, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 984, 504, 684, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1344, cbuffer, 496, 586, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1614, cbuffer, 586, 736, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2064, 1344, 1614, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2604, cbuffer, 961, 997, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2712, cbuffer, 997, 1057, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2892, 2604, 2712, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3108, cbuffer, 1147, 1207, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3288, cbuffer, 1207, 1307, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3588, 3108, 3288, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3948, cbuffer, 1457, 1547, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4218, cbuffer, 1547, 1697, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4668, 3948, 4218, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5208, cbuffer, 1922, 2048, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5586, cbuffer, 2048, 2258, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 6216, 5208, 5586, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 6972, cbuffer, 2573, 2741, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7476, cbuffer, 2741, 3021, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8316, 6972, 7476, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 288, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 150, ckbuffer, 984, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 400, ckbuffer, 2064, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 2875, ckbuffer, 2892, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3025, ckbuffer, 3588, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3275, ckbuffer, 4668, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3650, ckbuffer, 6216, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 4175, ckbuffer, 8316, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 775, 0, 150, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 1225, 150, 400, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 1975, 775, 1225, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 4875, 2875, 3025, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5325, 3025, 3275, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 6075, 3275, 3650, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 7200, 3650, 4175, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 8775, 0, 4875, 5325, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 10125, 150, 5325, 6075, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 12375, 400, 6075, 7200, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 15750, 775, 8775, 10125, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 18450, 1225, 10125, 12375, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_fdxx(skbuffer, 22950, 1975, 15750, 18450, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 22950, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 875, skbuffer, 24450, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1750, skbuffer, 25950, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1000RecFDDD_hpp */
