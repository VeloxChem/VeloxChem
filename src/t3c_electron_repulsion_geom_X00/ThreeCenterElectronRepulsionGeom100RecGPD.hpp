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

#ifndef ThreeCenterElectronRepulsionGeom100RecGPD_hpp
#define ThreeCenterElectronRepulsionGeom100RecGPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecGXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
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

/// @brief Computes d^(1)/dA^(1)(G|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_gpd(T& distributor,
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

    CSimdArray<double> pbuffer(1880, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1216, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(918, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(405, 1);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 115, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 118, 3, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 121, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 124, 2, 14, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 133, 3, 17, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 142, 4, 20, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 151, 8, 29, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 169, 11, 35, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 187, 14, 41, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 205, 17, 47, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 223, 20, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 241, 29, 65, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 271, 35, 75, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 301, 41, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 331, 47, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 361, 53, 105, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 391, 2, 3, 121, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 397, 8, 11, 115, 124, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 415, 11, 14, 118, 133, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 433, 14, 17, 121, 142, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 451, 29, 35, 124, 187, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 487, 35, 41, 133, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 523, 41, 47, 142, 223, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 559, 65, 75, 187, 301, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 619, 75, 85, 205, 331, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 679, 85, 95, 223, 361, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 739, 115, 118, 391, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 749, 124, 133, 391, 433, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 779, 151, 169, 397, 451, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 839, 169, 187, 415, 487, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 899, 187, 205, 433, 523, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 959, 241, 271, 451, 559, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1059, 271, 301, 487, 619, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1159, 301, 331, 523, 679, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 1259, 397, 415, 739, 749, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 1304, 451, 487, 749, 899, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 1394, 559, 619, 899, 1159, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsd(pbuffer, 1544, 779, 839, 1259, 1304, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 1670, 959, 1059, 1304, 1394, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 779, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 959, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1544, 1670});

                pbuffer.scale(2.0 * a_exp, {1670, 1880});

                t2cfunc::reduce(cbuffer, 880, pbuffer, 1544, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1006, pbuffer, 1670, 210, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 160, 0, 880, 0, 2);

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 430, 60, 1006, 0, 3);

            t3cfunc::bra_transform<4>(skbuffer, 0, cbuffer, 160, 0, 2);

            t3cfunc::bra_transform<4>(skbuffer, 54, cbuffer, 250, 0, 2);

            t3cfunc::bra_transform<4>(skbuffer, 108, cbuffer, 340, 0, 2);

            t3cfunc::bra_transform<4>(skbuffer, 162, cbuffer, 430, 0, 3);

            t3cfunc::bra_transform<4>(skbuffer, 252, cbuffer, 580, 0, 3);

            t3cfunc::bra_transform<4>(skbuffer, 342, cbuffer, 730, 0, 3);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 432, 0, 162, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 594, 54, 252, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 756, 108, 342, cfactors, 6, 4);

            t3cfunc::ket_transform<1, 2>(sbuffer, 0, skbuffer, 432, 4);

            t3cfunc::ket_transform<1, 2>(sbuffer, 135, skbuffer, 594, 4);

            t3cfunc::ket_transform<1, 2>(sbuffer, 270, skbuffer, 756, 4);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 4, 1, 2, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecGPD_hpp */
