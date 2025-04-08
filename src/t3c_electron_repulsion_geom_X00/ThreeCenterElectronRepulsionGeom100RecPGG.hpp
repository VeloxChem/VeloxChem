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

#ifndef ThreeCenterElectronRepulsionGeom100RecPGG_hpp
#define ThreeCenterElectronRepulsionGeom100RecPGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecPXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)(P|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_pgg(T& distributor,
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

    CSimdArray<double> pfactors(29, ket_npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(2160, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2320, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(12726, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(729, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

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

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_a);

                t3cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 11);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 8, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 9, pfactors, 16, bf_data, 9);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 10, pfactors, 16, bf_data, 10);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 280, 95, 105, 175, 190, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 301, 105, 115, 190, 205, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 322, 115, 125, 205, 220, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 343, 125, 135, 220, 235, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 364, 135, 145, 235, 250, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 385, 145, 155, 250, 265, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 406, 175, 190, 280, 301, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 434, 190, 205, 301, 322, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 462, 205, 220, 322, 343, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 490, 220, 235, 343, 364, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 518, 235, 250, 364, 385, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 546, 280, 301, 406, 434, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 582, 301, 322, 434, 462, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 618, 322, 343, 462, 490, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 654, 343, 364, 490, 518, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 690, 406, 434, 546, 582, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 735, 434, 462, 582, 618, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 780, 462, 490, 618, 654, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 825, 53, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 855, 115, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 900, 205, 322, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 963, 322, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1047, 462, 618, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1155, 618, 780, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1290, 175, 190, 825, 855, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1380, 280, 301, 855, 900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 1506, 406, 434, 900, 963, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 1674, 546, 582, 963, 1047, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 1890, 690, 735, 1047, 1155, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 546, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 690, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1290, 1380});

                pbuffer.scale(2.0 * a_exp, {1380, 1506});

                pbuffer.scale(2.0 * a_exp, {1506, 1674});

                pbuffer.scale(2.0 * a_exp, {1674, 1890});

                pbuffer.scale(2.0 * a_exp, {1890, 2160});

                t2cfunc::reduce(cbuffer, 1450, pbuffer, 1290, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1540, pbuffer, 1380, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1666, pbuffer, 1506, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1834, pbuffer, 1674, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2050, pbuffer, 1890, 270, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 145, 0, 1450, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 280, 15, 1540, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 469, 36, 1666, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 721, 64, 1834, 0, 7);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 1045, 100, 2050, 0, 8);

            t3cfunc::bra_transform<1>(skbuffer, 0, cbuffer, 145, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 45, cbuffer, 190, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 90, cbuffer, 235, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 135, cbuffer, 280, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 198, cbuffer, 343, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 261, cbuffer, 406, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 324, cbuffer, 469, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 408, cbuffer, 553, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 492, cbuffer, 637, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 576, cbuffer, 721, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 684, cbuffer, 829, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 792, cbuffer, 937, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 900, cbuffer, 1045, 0, 8);

            t3cfunc::bra_transform<1>(skbuffer, 1035, cbuffer, 1180, 0, 8);

            t3cfunc::bra_transform<1>(skbuffer, 1170, cbuffer, 1315, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1305, 0, 135, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1440, 45, 198, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1575, 90, 261, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1710, 135, 324, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1899, 198, 408, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2088, 261, 492, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2277, 324, 576, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2529, 408, 684, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2781, 492, 792, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 3033, 576, 900, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 3357, 684, 1035, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 3681, 792, 1170, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4005, 1305, 1710, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4275, 1440, 1899, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4545, 1575, 2088, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 4815, 1710, 2277, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 5193, 1899, 2529, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 5571, 2088, 2781, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 5949, 2277, 3033, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 6453, 2529, 3357, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 6957, 2781, 3681, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 7461, 4005, 4815, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 7911, 4275, 5193, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 8361, 4545, 5571, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 8811, 4815, 5949, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 9441, 5193, 6453, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 10071, 5571, 6957, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 10701, 7461, 8811, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 11376, 7911, 9441, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 12051, 8361, 10071, cfactors, 6, 1);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 10701, 1);

            t3cfunc::ket_transform<4, 4>(sbuffer, 243, skbuffer, 11376, 1);

            t3cfunc::ket_transform<4, 4>(sbuffer, 486, skbuffer, 12051, 1);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 1, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecPGG_hpp */
