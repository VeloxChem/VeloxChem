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

#ifndef ThreeCenterElectronRepulsionGeom010RecFGS_hpp
#define ThreeCenterElectronRepulsionGeom010RecFGS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
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

/// @brief Computes d^(1)/dC^(1)(F|1/|r-r'||GS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_fgs(T& distributor,
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

    CSimdArray<double> pbuffer(1638, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(760, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5390, 1);

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

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 115, 29, 35, 65, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 130, 35, 41, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 145, 41, 47, 85, 95, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 160, 47, 53, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 175, 65, 75, 115, 130, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 196, 75, 85, 130, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 217, 85, 95, 145, 160, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 238, 0, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 241, 1, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 244, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 247, 0, 8, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 256, 1, 11, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 265, 2, 14, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 274, 8, 29, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 292, 11, 35, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 310, 14, 41, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 328, 29, 65, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 358, 35, 75, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 388, 41, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 418, 65, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 463, 75, 130, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 508, 85, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 553, 115, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 616, 130, 196, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 679, 145, 217, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 742, 0, 1, 244, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 748, 8, 11, 244, 265, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 766, 29, 35, 265, 310, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 802, 65, 75, 310, 388, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 862, 115, 130, 388, 508, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 952, 175, 196, 508, 679, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 1078, 238, 241, 742, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 1088, 247, 256, 742, 748, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1118, 274, 292, 748, 766, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1178, 328, 358, 766, 802, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 1278, 418, 463, 802, 862, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 1428, 553, 616, 862, 952, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1078, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 1088, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 1118, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1178, 100, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {1078, 1088});

                pbuffer.scale(pfactors, 0, 2.0, {1088, 1118});

                pbuffer.scale(pfactors, 0, 2.0, {1118, 1178});

                pbuffer.scale(pfactors, 0, 2.0, {1178, 1278});

                pbuffer.scale(pfactors, 0, 2.0, {1278, 1428});

                pbuffer.scale(pfactors, 0, 2.0, {1428, 1638});

                t2cfunc::reduce(cbuffer, 200, pbuffer, 1078, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 210, pbuffer, 1088, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 240, pbuffer, 1118, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 1178, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 1278, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 1428, 210, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<3>(skbuffer, 0, cbuffer, 0, 0, 0);

            t3cfunc::bra_transform<3>(skbuffer, 7, cbuffer, 10, 0, 1);

            t3cfunc::bra_transform<3>(skbuffer, 28, cbuffer, 40, 0, 2);

            t3cfunc::bra_transform<3>(skbuffer, 70, cbuffer, 100, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 4263, cbuffer, 200, 0, 0);

            t3cfunc::bra_transform<3>(skbuffer, 4270, cbuffer, 210, 0, 1);

            t3cfunc::bra_transform<3>(skbuffer, 4291, cbuffer, 240, 0, 2);

            t3cfunc::bra_transform<3>(skbuffer, 4333, cbuffer, 300, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 4403, cbuffer, 400, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 4508, cbuffer, 550, 0, 5);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 140, 0, 7, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 224, 7, 28, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 476, 28, 70, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xds(skbuffer, 1610, 140, 224, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdp(skbuffer, 1778, 224, 476, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xfs(skbuffer, 3038, 1610, 1778, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xss(skbuffer, 4655, 4263, 4270, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 4676, 4270, 4291, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 4739, 4291, 4333, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 4865, 4333, 4403, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 5075, 4403, 4508, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xps(skbuffer, 161, 0, 4655, 4676, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 287, 7, 4676, 4739, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 602, 28, 4739, 4865, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 980, 70, 4865, 5075, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xds(skbuffer, 1652, 140, 161, 287, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 1904, 224, 287, 602, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 2282, 476, 602, 980, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfs(skbuffer, 3108, 1610, 1652, 1904, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfp(skbuffer, 3318, 1778, 1904, 2282, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xgs(skbuffer, 3948, 3038, 3108, 3318, cfactors, 6, 3);

            t3cfunc::ket_transform<4, 0>(sbuffer, 0, skbuffer, 3948, 3);

            t3cfunc::ket_transform<4, 0>(sbuffer, 63, skbuffer, 4053, 3);

            t3cfunc::ket_transform<4, 0>(sbuffer, 126, skbuffer, 4158, 3);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 3, 4, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecFGS_hpp */
