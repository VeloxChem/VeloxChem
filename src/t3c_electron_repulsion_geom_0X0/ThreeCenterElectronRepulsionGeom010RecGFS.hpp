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

#ifndef ThreeCenterElectronRepulsionGeom010RecGFS_hpp
#define ThreeCenterElectronRepulsionGeom010RecGFS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
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

/// @brief Computes d^(1)/dC^(1)(G|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_gfs(T& distributor,
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

    CSimdArray<double> pbuffer(2030, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(675, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2835, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(189, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 9);

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

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 210, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 213, 3, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 216, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 219, 2, 15, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 228, 3, 18, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 237, 4, 21, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 246, 15, 45, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 264, 18, 51, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 282, 21, 57, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 300, 45, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 330, 51, 105, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 360, 57, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 390, 95, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 435, 105, 180, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 480, 115, 195, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 525, 0, 1, 210, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 531, 1, 2, 213, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 537, 2, 3, 216, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 543, 9, 12, 210, 219, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 561, 12, 15, 213, 228, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 579, 15, 18, 216, 237, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 597, 33, 39, 219, 246, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 633, 39, 45, 228, 264, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 669, 45, 51, 237, 282, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 705, 75, 85, 246, 300, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 765, 85, 95, 264, 330, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 825, 95, 105, 282, 360, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 885, 135, 150, 300, 390, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 975, 150, 165, 330, 435, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1065, 165, 180, 360, 480, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 1155, 210, 213, 537, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 1165, 219, 228, 537, 579, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1195, 246, 264, 579, 669, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1255, 300, 330, 669, 825, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 1355, 390, 435, 825, 1065, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 1505, 525, 531, 1155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 1520, 543, 561, 1155, 1165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 1565, 597, 633, 1165, 1195, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 1655, 705, 765, 1195, 1255, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 1805, 885, 975, 1255, 1355, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1505, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 1520, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 1565, 90, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {1505, 1520});

                pbuffer.scale(pfactors, 0, 2.0, {1520, 1565});

                pbuffer.scale(pfactors, 0, 2.0, {1565, 1655});

                pbuffer.scale(pfactors, 0, 2.0, {1655, 1805});

                pbuffer.scale(pfactors, 0, 2.0, {1805, 2030});

                t2cfunc::reduce(cbuffer, 150, pbuffer, 1505, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 165, pbuffer, 1520, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 210, pbuffer, 1565, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 1655, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 1805, 225, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<4>(skbuffer, 0, cbuffer, 0, 0, 0);

            t3cfunc::bra_transform<4>(skbuffer, 9, cbuffer, 15, 0, 1);

            t3cfunc::bra_transform<4>(skbuffer, 36, cbuffer, 60, 0, 2);

            t3cfunc::bra_transform<4>(skbuffer, 1980, cbuffer, 150, 0, 0);

            t3cfunc::bra_transform<4>(skbuffer, 1989, cbuffer, 165, 0, 1);

            t3cfunc::bra_transform<4>(skbuffer, 2016, cbuffer, 210, 0, 2);

            t3cfunc::bra_transform<4>(skbuffer, 2070, cbuffer, 300, 0, 3);

            t3cfunc::bra_transform<4>(skbuffer, 2160, cbuffer, 450, 0, 4);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 90, 0, 9, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 198, 9, 36, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xds(skbuffer, 1008, 90, 198, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xss(skbuffer, 2295, 1980, 1989, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 2322, 1989, 2016, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 2403, 2016, 2070, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 2565, 2070, 2160, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xps(skbuffer, 117, 0, 2295, 2322, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 279, 9, 2322, 2403, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 522, 36, 2403, 2565, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xds(skbuffer, 1062, 90, 117, 279, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 1224, 198, 279, 522, cfactors, 6, 4);

            t3ceri::comp_ket_geom010_electron_repulsion_xfs(skbuffer, 1710, 1008, 1062, 1224, cfactors, 6, 4);

            t3cfunc::ket_transform<3, 0>(sbuffer, 0, skbuffer, 1710, 4);

            t3cfunc::ket_transform<3, 0>(sbuffer, 63, skbuffer, 1800, 4);

            t3cfunc::ket_transform<3, 0>(sbuffer, 126, skbuffer, 1890, 4);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 4, 3, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecGFS_hpp */
