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

#ifndef ThreeCenterElectronRepulsionGeom010RecPGP_hpp
#define ThreeCenterElectronRepulsionGeom010RecPGP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
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

/// @brief Computes d^(1)/dC^(1)(P|1/|r-r'||GP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_pgp(T& distributor,
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

    CSimdArray<double> pbuffer(459, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(351, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4428, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(243, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 8);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 95, 25, 31, 55, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 110, 31, 37, 65, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 125, 37, 43, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 140, 55, 65, 95, 110, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 161, 65, 75, 110, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 182, 95, 110, 140, 161, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 210, 0, 7, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 219, 7, 25, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 237, 25, 55, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 267, 55, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 312, 95, 140, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 375, 140, 182, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 210, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 219, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 237, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 57, pbuffer, 267, 45, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {210, 219});

                pbuffer.scale(pfactors, 0, 2.0, {219, 237});

                pbuffer.scale(pfactors, 0, 2.0, {237, 267});

                pbuffer.scale(pfactors, 0, 2.0, {267, 312});

                pbuffer.scale(pfactors, 0, 2.0, {312, 375});

                pbuffer.scale(pfactors, 0, 2.0, {375, 459});

                t2cfunc::reduce(cbuffer, 102, pbuffer, 210, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 219, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 129, pbuffer, 237, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 159, pbuffer, 267, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 312, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 267, pbuffer, 375, 84, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<1>(skbuffer, 0, cbuffer, 0, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 9, cbuffer, 9, 0, 2);

            t3cfunc::bra_transform<1>(skbuffer, 27, cbuffer, 27, 0, 3);

            t3cfunc::bra_transform<1>(skbuffer, 57, cbuffer, 57, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 3684, cbuffer, 102, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 3693, cbuffer, 111, 0, 2);

            t3cfunc::bra_transform<1>(skbuffer, 3711, cbuffer, 129, 0, 3);

            t3cfunc::bra_transform<1>(skbuffer, 3741, cbuffer, 159, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 3786, cbuffer, 204, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 3849, cbuffer, 267, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 102, 0, 9, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 210, 9, 27, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 426, 27, 57, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdp(skbuffer, 1191, 102, 210, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdd(skbuffer, 1407, 210, 426, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfp(skbuffer, 2379, 1191, 1407, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 3933, 3684, 3693, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 3960, 3693, 3711, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 4014, 3711, 3741, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 4104, 3741, 3786, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsh(skbuffer, 4239, 3786, 3849, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 129, 0, 3933, 3960, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 264, 9, 3960, 4014, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 516, 27, 4014, 4104, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 786, 57, 4104, 4239, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 1245, 102, 129, 264, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 1515, 210, 264, 516, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xdf(skbuffer, 1839, 426, 516, 786, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xfp(skbuffer, 2469, 1191, 1245, 1515, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xfd(skbuffer, 2739, 1407, 1515, 1839, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xgp(skbuffer, 3279, 2379, 2469, 2739, cfactors, 6, 1);

            t3cfunc::ket_transform<4, 1>(sbuffer, 0, skbuffer, 3279, 1);

            t3cfunc::ket_transform<4, 1>(sbuffer, 81, skbuffer, 3414, 1);

            t3cfunc::ket_transform<4, 1>(sbuffer, 162, skbuffer, 3549, 1);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 1, 4, 1, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecPGP_hpp */
