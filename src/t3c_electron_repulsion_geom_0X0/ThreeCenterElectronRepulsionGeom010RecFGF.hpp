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

#ifndef ThreeCenterElectronRepulsionGeom010RecFGF_hpp
#define ThreeCenterElectronRepulsionGeom010RecFGF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
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

/// @brief Computes d^(1)/dC^(1)(F|1/|r-r'||GF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_fgf(T& distributor,
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

    CSimdArray<double> pbuffer(4763, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2290, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(25487, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

    // setup Boys fuction data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 12);

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

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 825, 2, 17, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 834, 17, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 852, 41, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 882, 47, 105, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 912, 53, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 942, 95, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 987, 105, 190, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1032, 115, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1077, 175, 280, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1140, 190, 301, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1203, 205, 322, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1266, 280, 406, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1350, 301, 434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1434, 322, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1518, 406, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1626, 434, 582, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1734, 462, 618, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1842, 546, 690, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1977, 582, 735, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2112, 618, 780, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2247, 41, 47, 825, 834, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2283, 95, 105, 834, 912, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2343, 175, 190, 912, 1032, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2433, 280, 301, 1032, 1203, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2559, 406, 434, 1203, 1434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 2727, 546, 582, 1434, 1734, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 2943, 690, 735, 1734, 2112, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 3213, 852, 882, 2247, 2283, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 3313, 942, 987, 2283, 2343, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 3463, 1077, 1140, 2343, 2433, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 3673, 1266, 1350, 2433, 2559, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 3953, 1518, 1626, 2559, 2727, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 4313, 1842, 1977, 2727, 2943, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3213, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 3313, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 3463, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 3673, 280, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {3213, 3313});

                pbuffer.scale(pfactors, 0, 2.0, {3313, 3463});

                pbuffer.scale(pfactors, 0, 2.0, {3463, 3673});

                pbuffer.scale(pfactors, 0, 2.0, {3673, 3953});

                pbuffer.scale(pfactors, 0, 2.0, {3953, 4313});

                pbuffer.scale(pfactors, 0, 2.0, {4313, 4763});

                t2cfunc::reduce(cbuffer, 740, pbuffer, 3213, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 3313, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 990, pbuffer, 3463, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1200, pbuffer, 3673, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1480, pbuffer, 3953, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1840, pbuffer, 4313, 450, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<3>(skbuffer, 0, cbuffer, 0, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 70, cbuffer, 100, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 175, cbuffer, 250, 0, 5);

            t3cfunc::bra_transform<3>(skbuffer, 322, cbuffer, 460, 0, 6);

            t3cfunc::bra_transform<3>(skbuffer, 22092, cbuffer, 740, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 22162, cbuffer, 840, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 22267, cbuffer, 990, 0, 5);

            t3cfunc::bra_transform<3>(skbuffer, 22414, cbuffer, 1200, 0, 6);

            t3cfunc::bra_transform<3>(skbuffer, 22610, cbuffer, 1480, 0, 7);

            t3cfunc::bra_transform<3>(skbuffer, 22862, cbuffer, 1840, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 518, 0, 70, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1358, 70, 175, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2618, 175, 322, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 6146, 518, 1358, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 7826, 1358, 2618, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 12992, 6146, 7826, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 23177, 22092, 22162, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 23387, 22162, 22267, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsh(skbuffer, 23702, 22267, 22414, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsi(skbuffer, 24143, 22414, 22610, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsk(skbuffer, 24731, 22610, 22862, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 728, 0, 23177, 23387, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 1673, 70, 23387, 23702, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 3059, 175, 23702, 24143, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpi(skbuffer, 4382, 322, 24143, 24731, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdf(skbuffer, 6566, 518, 728, 1673, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 8456, 1358, 1673, 3059, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdh(skbuffer, 10346, 2618, 3059, 4382, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xff(skbuffer, 13692, 6146, 6566, 8456, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfg(skbuffer, 15792, 7826, 8456, 10346, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xgf(skbuffer, 18942, 12992, 13692, 15792, cfactors, 6, 3);

            t3cfunc::ket_transform<4, 3>(sbuffer, 0, skbuffer, 18942, 3);

            t3cfunc::ket_transform<4, 3>(sbuffer, 441, skbuffer, 19992, 3);

            t3cfunc::ket_transform<4, 3>(sbuffer, 882, skbuffer, 21042, 3);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 3, 4, 3, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecFGF_hpp */
