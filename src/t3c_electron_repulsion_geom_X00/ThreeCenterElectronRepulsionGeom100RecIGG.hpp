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

#ifndef ThreeCenterElectronRepulsionGeom100RecIGG_hpp
#define ThreeCenterElectronRepulsionGeom100RecIGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDI.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXGG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPK.hpp"
#include "ThreeCenterElectronRepulsionGeom100ContrRecIXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)(I|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_igg(T& distributor,
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

    CSimdArray<double> pbuffer(43535, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(25665, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(55146, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3159, 1);

    // setup Boys fuction data

    const CBoysFunc<15> bf_table;

    CSimdArray<double> bf_data(17, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 16, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 16);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 8, pfactors, 16, bf_data, 9);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 9, pfactors, 16, bf_data, 10);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 10, pfactors, 16, bf_data, 11);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 11, pfactors, 16, bf_data, 12);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 12, pfactors, 16, bf_data, 13);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 13, pfactors, 16, bf_data, 14);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 14, pfactors, 16, bf_data, 15);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 33, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 36, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 39, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 42, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 45, 10, 11, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 48, 11, 12, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 51, 12, 13, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 54, 13, 14, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 0, 1, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 1, 2, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 2, 3, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 75, 3, 4, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 81, 4, 5, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 87, 5, 6, 30, 33, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 93, 6, 7, 33, 36, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 99, 7, 8, 36, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 105, 8, 9, 39, 42, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 111, 9, 10, 42, 45, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 117, 10, 11, 45, 48, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 123, 11, 12, 48, 51, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 129, 12, 13, 51, 54, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 15, 18, 57, 63, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 18, 21, 63, 69, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 21, 24, 69, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 24, 27, 75, 81, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 175, 27, 30, 81, 87, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 185, 30, 33, 87, 93, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 195, 33, 36, 93, 99, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 205, 36, 39, 99, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 215, 39, 42, 105, 111, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 225, 42, 45, 111, 117, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 235, 45, 48, 117, 123, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 245, 48, 51, 123, 129, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 255, 57, 63, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 270, 63, 69, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 285, 69, 75, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 300, 75, 81, 165, 175, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 315, 81, 87, 175, 185, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 330, 87, 93, 185, 195, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 345, 93, 99, 195, 205, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 360, 99, 105, 205, 215, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 375, 105, 111, 215, 225, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 390, 111, 117, 225, 235, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 405, 117, 123, 235, 245, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 420, 135, 145, 255, 270, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 441, 145, 155, 270, 285, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 462, 155, 165, 285, 300, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 483, 165, 175, 300, 315, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 504, 175, 185, 315, 330, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 525, 185, 195, 330, 345, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 546, 195, 205, 345, 360, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 567, 205, 215, 360, 375, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 588, 215, 225, 375, 390, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 609, 225, 235, 390, 405, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 630, 255, 270, 420, 441, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 658, 270, 285, 441, 462, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 686, 285, 300, 462, 483, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 714, 300, 315, 483, 504, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 742, 315, 330, 504, 525, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 770, 330, 345, 525, 546, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 798, 345, 360, 546, 567, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 826, 360, 375, 567, 588, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 854, 375, 390, 588, 609, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 882, 420, 441, 630, 658, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 918, 441, 462, 658, 686, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 954, 462, 483, 686, 714, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 990, 483, 504, 714, 742, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1026, 504, 525, 742, 770, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1062, 525, 546, 770, 798, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1098, 546, 567, 798, 826, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1134, 567, 588, 826, 854, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1170, 630, 658, 882, 918, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1215, 658, 686, 918, 954, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1260, 686, 714, 954, 990, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1305, 714, 742, 990, 1026, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1350, 742, 770, 1026, 1062, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1395, 770, 798, 1062, 1098, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1440, 798, 826, 1098, 1134, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 1485, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 1488, 5, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 1491, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1494, 4, 27, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1503, 5, 30, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1512, 6, 33, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1521, 21, 69, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1539, 24, 75, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1557, 27, 81, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1575, 30, 87, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1593, 33, 93, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1611, 69, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1641, 75, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1671, 81, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1701, 87, 185, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1731, 93, 195, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1761, 135, 255, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1806, 145, 270, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1851, 155, 285, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1896, 165, 300, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1941, 175, 315, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1986, 185, 330, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 2031, 195, 345, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2076, 255, 420, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2139, 270, 441, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2202, 285, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2265, 300, 483, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2328, 315, 504, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2391, 330, 525, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2454, 345, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2517, 420, 630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2601, 441, 658, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2685, 462, 686, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2769, 483, 714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2853, 504, 742, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2937, 525, 770, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 3021, 546, 798, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3105, 630, 882, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3213, 658, 918, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3321, 686, 954, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3429, 714, 990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3537, 742, 1026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3645, 770, 1062, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3753, 798, 1098, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3861, 882, 1170, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3996, 918, 1215, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 4131, 954, 1260, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 4266, 990, 1305, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 4401, 1026, 1350, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 4536, 1062, 1395, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 4671, 1098, 1440, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 4806, 4, 5, 1491, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 4812, 21, 24, 1485, 1494, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 4830, 24, 27, 1488, 1503, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 4848, 27, 30, 1491, 1512, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 4866, 69, 75, 1494, 1557, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 4902, 75, 81, 1503, 1575, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 4938, 81, 87, 1512, 1593, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 4974, 135, 145, 1521, 1611, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 5034, 145, 155, 1539, 1641, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 5094, 155, 165, 1557, 1671, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 5154, 165, 175, 1575, 1701, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 5214, 175, 185, 1593, 1731, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 5274, 255, 270, 1611, 1851, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 5364, 270, 285, 1641, 1896, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 5454, 285, 300, 1671, 1941, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 5544, 300, 315, 1701, 1986, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 5634, 315, 330, 1731, 2031, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 5724, 420, 441, 1851, 2202, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 5850, 441, 462, 1896, 2265, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 5976, 462, 483, 1941, 2328, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 6102, 483, 504, 1986, 2391, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 6228, 504, 525, 2031, 2454, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 6354, 630, 658, 2202, 2685, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 6522, 658, 686, 2265, 2769, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 6690, 686, 714, 2328, 2853, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 6858, 714, 742, 2391, 2937, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 7026, 742, 770, 2454, 3021, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 7194, 882, 918, 2685, 3321, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 7410, 918, 954, 2769, 3429, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 7626, 954, 990, 2853, 3537, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 7842, 990, 1026, 2937, 3645, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 8058, 1026, 1062, 3021, 3753, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 8274, 1170, 1215, 3321, 4131, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 8544, 1215, 1260, 3429, 4266, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 8814, 1260, 1305, 3537, 4401, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 9084, 1305, 1350, 3645, 4536, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 9354, 1350, 1395, 3753, 4671, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 9624, 1485, 1488, 4806, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 9634, 1494, 1503, 4806, 4848, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 9664, 1521, 1539, 4812, 4866, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 9724, 1539, 1557, 4830, 4902, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 9784, 1557, 1575, 4848, 4938, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 9844, 1611, 1641, 4866, 5094, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 9944, 1641, 1671, 4902, 5154, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 10044, 1671, 1701, 4938, 5214, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 10144, 1761, 1806, 4974, 5274, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 10294, 1806, 1851, 5034, 5364, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 10444, 1851, 1896, 5094, 5454, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 10594, 1896, 1941, 5154, 5544, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 10744, 1941, 1986, 5214, 5634, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 10894, 2076, 2139, 5274, 5724, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 11104, 2139, 2202, 5364, 5850, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 11314, 2202, 2265, 5454, 5976, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 11524, 2265, 2328, 5544, 6102, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 11734, 2328, 2391, 5634, 6228, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 11944, 2517, 2601, 5724, 6354, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 12224, 2601, 2685, 5850, 6522, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 12504, 2685, 2769, 5976, 6690, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 12784, 2769, 2853, 6102, 6858, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 13064, 2853, 2937, 6228, 7026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 13344, 3105, 3213, 6354, 7194, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 13704, 3213, 3321, 6522, 7410, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 14064, 3321, 3429, 6690, 7626, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 14424, 3429, 3537, 6858, 7842, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 14784, 3537, 3645, 7026, 8058, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 15144, 3861, 3996, 7194, 8274, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 15594, 3996, 4131, 7410, 8544, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 16044, 4131, 4266, 7626, 8814, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 16494, 4266, 4401, 7842, 9084, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 16944, 4401, 4536, 8058, 9354, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 17394, 4812, 4830, 9624, 9634, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 17439, 4866, 4902, 9634, 9784, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 17529, 4974, 5034, 9664, 9844, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 17679, 5034, 5094, 9724, 9944, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 17829, 5094, 5154, 9784, 10044, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 17979, 5274, 5364, 9844, 10444, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 18204, 5364, 5454, 9944, 10594, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 18429, 5454, 5544, 10044, 10744, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 18654, 5724, 5850, 10444, 11314, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 18969, 5850, 5976, 10594, 11524, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 19284, 5976, 6102, 10744, 11734, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 19599, 6354, 6522, 11314, 12504, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 20019, 6522, 6690, 11524, 12784, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 20439, 6690, 6858, 11734, 13064, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 20859, 7194, 7410, 12504, 14064, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 21399, 7410, 7626, 12784, 14424, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 21939, 7626, 7842, 13064, 14784, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 22479, 8274, 8544, 14064, 16044, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 23154, 8544, 8814, 14424, 16494, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 23829, 8814, 9084, 14784, 16944, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsd(pbuffer, 24504, 9664, 9724, 17394, 17439, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 24630, 9844, 9944, 17439, 17829, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 24840, 10144, 10294, 17529, 17979, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 25155, 10294, 10444, 17679, 18204, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 25470, 10444, 10594, 17829, 18429, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 25785, 10894, 11104, 17979, 18654, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 26226, 11104, 11314, 18204, 18969, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 26667, 11314, 11524, 18429, 19284, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 27108, 11944, 12224, 18654, 19599, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 27696, 12224, 12504, 18969, 20019, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 28284, 12504, 12784, 19284, 20439, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsk(pbuffer, 28872, 13344, 13704, 19599, 20859, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsk(pbuffer, 29628, 13704, 14064, 20019, 21399, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsk(pbuffer, 30384, 14064, 14424, 20439, 21939, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsl(pbuffer, 31140, 15144, 15594, 20859, 22479, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsl(pbuffer, 32085, 15594, 16044, 21399, 23154, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsl(pbuffer, 33030, 16044, 16494, 21939, 23829, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isf(pbuffer, 33975, 17529, 17679, 24504, 24630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isg(pbuffer, 34255, 17979, 18204, 24630, 25470, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ish(pbuffer, 34675, 18654, 18969, 25470, 26667, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isi(pbuffer, 35263, 19599, 20019, 26667, 28284, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isk(pbuffer, 36047, 20859, 21399, 28284, 30384, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isl(pbuffer, 37055, 22479, 23154, 30384, 33030, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksg(pbuffer, 38315, 24840, 25155, 33975, 34255, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksh(pbuffer, 38855, 25785, 26226, 34255, 34675, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksi(pbuffer, 39611, 27108, 27696, 34675, 35263, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksk(pbuffer, 40619, 28872, 29628, 35263, 36047, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksl(pbuffer, 41915, 31140, 32085, 36047, 37055, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 24840, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 315, pbuffer, 25785, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 756, pbuffer, 27108, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1344, pbuffer, 28872, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2100, pbuffer, 31140, 945, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15225, pbuffer, 38315, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15765, pbuffer, 38855, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16521, pbuffer, 39611, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17529, pbuffer, 40619, 1296, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18825, pbuffer, 41915, 1620, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {38315, 38855});

                pbuffer.scale(2.0 * a_exp, {38855, 39611});

                pbuffer.scale(2.0 * a_exp, {39611, 40619});

                pbuffer.scale(2.0 * a_exp, {40619, 41915});

                pbuffer.scale(2.0 * a_exp, {41915, 43535});

                t2cfunc::reduce(cbuffer, 20445, pbuffer, 38315, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20985, pbuffer, 38855, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21741, pbuffer, 39611, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22749, pbuffer, 40619, 1296, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24045, pbuffer, 41915, 1620, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 3045, 0, 20445, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 4305, 315, 20985, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 6069, 756, 21741, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 8421, 1344, 22749, 0, 7);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 11445, 2100, 24045, 0, 8);

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 3045, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 195, cbuffer, 3465, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 390, cbuffer, 3885, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 585, cbuffer, 4305, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 858, cbuffer, 4893, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 1131, cbuffer, 5481, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 1404, cbuffer, 6069, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 1768, cbuffer, 6853, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 2132, cbuffer, 7637, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 2496, cbuffer, 8421, 0, 7);

            t3cfunc::bra_transform<6>(skbuffer, 2964, cbuffer, 9429, 0, 7);

            t3cfunc::bra_transform<6>(skbuffer, 3432, cbuffer, 10437, 0, 7);

            t3cfunc::bra_transform<6>(skbuffer, 3900, cbuffer, 11445, 0, 8);

            t3cfunc::bra_transform<6>(skbuffer, 4485, cbuffer, 12705, 0, 8);

            t3cfunc::bra_transform<6>(skbuffer, 5070, cbuffer, 13965, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 5655, 0, 585, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 6240, 195, 858, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 6825, 390, 1131, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 7410, 585, 1404, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 8229, 858, 1768, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 9048, 1131, 2132, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 9867, 1404, 2496, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 10959, 1768, 2964, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 12051, 2132, 3432, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 13143, 2496, 3900, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 14547, 2964, 4485, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 15951, 3432, 5070, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 17355, 5655, 7410, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 18525, 6240, 8229, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 19695, 6825, 9048, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 20865, 7410, 9867, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 22503, 8229, 10959, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 24141, 9048, 12051, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 25779, 9867, 13143, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 27963, 10959, 14547, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 30147, 12051, 15951, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 32331, 17355, 20865, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 34281, 18525, 22503, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 36231, 19695, 24141, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 38181, 20865, 25779, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 40911, 22503, 27963, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 43641, 24141, 30147, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 46371, 32331, 38181, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 49296, 34281, 40911, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 52221, 36231, 43641, cfactors, 6, 6);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 46371, 6);

            t3cfunc::ket_transform<4, 4>(sbuffer, 1053, skbuffer, 49296, 6);

            t3cfunc::ket_transform<4, 4>(sbuffer, 2106, skbuffer, 52221, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecIGG_hpp */
