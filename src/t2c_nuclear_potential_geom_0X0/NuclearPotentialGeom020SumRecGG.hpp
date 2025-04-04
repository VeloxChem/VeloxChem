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

#ifndef NuclearPotentialGeom020SumRecGG_hpp
#define NuclearPotentialGeom020SumRecGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDF.hpp"
#include "NuclearPotentialGeom010PrimRecDG.hpp"
#include "NuclearPotentialGeom010PrimRecFG.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPF.hpp"
#include "NuclearPotentialGeom010PrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSF.hpp"
#include "NuclearPotentialGeom010PrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialGeom020PrimRecDD.hpp"
#include "NuclearPotentialGeom020PrimRecDF.hpp"
#include "NuclearPotentialGeom020PrimRecDG.hpp"
#include "NuclearPotentialGeom020PrimRecFF.hpp"
#include "NuclearPotentialGeom020PrimRecFG.hpp"
#include "NuclearPotentialGeom020PrimRecGG.hpp"
#include "NuclearPotentialGeom020PrimRecPD.hpp"
#include "NuclearPotentialGeom020PrimRecPF.hpp"
#include "NuclearPotentialGeom020PrimRecPG.hpp"
#include "NuclearPotentialGeom020PrimRecPP.hpp"
#include "NuclearPotentialGeom020PrimRecSD.hpp"
#include "NuclearPotentialGeom020PrimRecSF.hpp"
#include "NuclearPotentialGeom020PrimRecSG.hpp"
#include "NuclearPotentialGeom020PrimRecSP.hpp"
#include "NuclearPotentialGeom020PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (G|AG(2)|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_gg(T&                               distributor,
                                       const CGtoBlock&                 bra_gto_block,
                                       const CGtoBlock&                 ket_gto_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices,
                                       const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external quadrupoles data

    const auto quadrupoles = distributor.data();

    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.coordinates();

    const auto ket_gto_exps = ket_gto_block.exponents();

    const auto ket_gto_norms = ket_gto_block.normalization_factors();

    const auto ket_gto_indices = ket_gto_block.orbital_indices();

    const auto ket_npgtos = ket_gto_block.number_of_primitives();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> factors(20, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(14502, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1350, 1);

    CSimdArray<double> sbuffer(486, 1);

    // setup Boys function data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        factors.load(ket_gto_exps, ket_range, 0, ket_npgtos);

        factors.load(ket_gto_norms, ket_range, 1, ket_npgtos);

        factors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        sbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        pbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over contracted basis functions on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            cbuffer.zero();

            sbuffer.zero();

            const auto r_a = bra_gto_coords[j];

            t2cfunc::comp_distances_ab(factors, 5, 2, r_a);

            for (size_t k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 11, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 11);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 9, 0, bf_data, 9, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 10, 0, bf_data, 10, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 11, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 14, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 17, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 20, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 23, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 26, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 29, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 32, 9, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 35, 1, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 41, 2, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 47, 3, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 53, 4, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 59, 5, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 65, 6, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 71, 7, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 77, 8, 9, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 83, 9, 10, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 89, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 92, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 95, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 98, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 101, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 104, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 107, 2, 11, 14, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 116, 3, 14, 17, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 125, 4, 17, 20, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 134, 5, 20, 23, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 143, 6, 23, 26, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 152, 7, 26, 29, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 161, 8, 29, 32, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 170, 11, 35, 41, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 188, 14, 41, 47, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 206, 17, 47, 53, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 224, 20, 53, 59, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 242, 23, 59, 65, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 260, 26, 65, 71, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 278, 29, 71, 77, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 296, 32, 77, 83, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 314, 2, 3, 89, 92, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 320, 3, 4, 92, 95, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 326, 4, 5, 95, 98, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 332, 5, 6, 98, 101, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 338, 6, 7, 101, 104, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 344, 11, 14, 89, 107, 116, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 362, 14, 17, 92, 116, 125, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 380, 17, 20, 95, 125, 134, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 398, 20, 23, 98, 134, 143, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 416, 23, 26, 101, 143, 152, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 434, 26, 29, 104, 152, 161, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 452, 35, 41, 107, 170, 188, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 488, 41, 47, 116, 188, 206, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 524, 47, 53, 125, 206, 224, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 560, 53, 59, 134, 224, 242, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 596, 59, 65, 143, 242, 260, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 632, 65, 71, 152, 260, 278, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 668, 71, 77, 161, 278, 296, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 704, 89, 92, 314, 320, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 714, 92, 95, 320, 326, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 724, 95, 98, 326, 332, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 734, 98, 101, 332, 338, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 744, 107, 116, 314, 344, 362, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 774, 116, 125, 320, 362, 380, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 804, 125, 134, 326, 380, 398, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 834, 134, 143, 332, 398, 416, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 864, 143, 152, 338, 416, 434, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 894, 170, 188, 344, 452, 488, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 954, 188, 206, 362, 488, 524, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 1014, 206, 224, 380, 524, 560, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 1074, 224, 242, 398, 560, 596, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 1134, 242, 260, 416, 596, 632, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 1194, 260, 278, 434, 632, 668, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 1254, 314, 320, 704, 714, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 1269, 320, 326, 714, 724, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 1284, 326, 332, 724, 734, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1299, 344, 362, 704, 744, 774, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1344, 362, 380, 714, 774, 804, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1389, 380, 398, 724, 804, 834, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1434, 398, 416, 734, 834, 864, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1479, 452, 488, 744, 894, 954, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1569, 488, 524, 774, 954, 1014, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1659, 524, 560, 804, 1014, 1074, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1749, 560, 596, 834, 1074, 1134, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1839, 596, 632, 864, 1134, 1194, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1929, 35, 41, 107, 170, 188, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1983, 41, 47, 116, 188, 206, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 2037, 47, 53, 125, 206, 224, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 2091, 53, 59, 134, 224, 242, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 2145, 107, 116, 314, 344, 362, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 2199, 116, 125, 320, 362, 380, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 2253, 125, 134, 326, 380, 398, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 2307, 170, 188, 344, 452, 488, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 2415, 188, 206, 362, 488, 524, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 2523, 206, 224, 380, 524, 560, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 2631, 224, 242, 398, 560, 596, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 2739, 314, 320, 704, 714, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 2769, 320, 326, 714, 724, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2799, 344, 362, 704, 744, 774, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2889, 362, 380, 714, 774, 804, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2979, 380, 398, 724, 804, 834, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 3069, 452, 488, 744, 894, 954, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 3249, 488, 524, 774, 954, 1014, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 3429, 524, 560, 804, 1014, 1074, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 3609, 560, 596, 834, 1074, 1134, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 3789, 704, 714, 1254, 1269, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 3834, 714, 724, 1269, 1284, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 3879, 744, 774, 1254, 1299, 1344, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 4014, 774, 804, 1269, 1344, 1389, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 4149, 804, 834, 1284, 1389, 1434, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 4284, 894, 954, 1299, 1479, 1569, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 4554, 954, 1014, 1344, 1569, 1659, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 4824, 1014, 1074, 1389, 1659, 1749, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 5094, 1074, 1134, 1434, 1749, 1839, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 5364, 452, 488, 1929, 1983, 2145, 2307, 2415, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 5580, 488, 524, 1983, 2037, 2199, 2415, 2523, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 5796, 524, 560, 2037, 2091, 2253, 2523, 2631, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 6012, 744, 774, 2145, 2199, 2739, 2799, 2889, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 6192, 774, 804, 2199, 2253, 2769, 2889, 2979, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 6372, 894, 954, 2307, 2415, 2799, 3069, 3249, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 6732, 954, 1014, 2415, 2523, 2889, 3249, 3429, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 7092, 1014, 1074, 2523, 2631, 2979, 3429, 3609, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 7452, 1254, 1269, 2739, 2769, 3789, 3834, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 7542, 1299, 1344, 2799, 2889, 3789, 3879, 4014, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 7812, 1344, 1389, 2889, 2979, 3834, 4014, 4149, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 8082, 1479, 1569, 3069, 3249, 3879, 4284, 4554, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 8622, 1569, 1659, 3249, 3429, 4014, 4554, 4824, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 9162, 1659, 1749, 3429, 3609, 4149, 4824, 5094, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ff(pbuffer, 9702, 3069, 3249, 5364, 5580, 6012, 6372, 6732, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ff(
                        pbuffer, 10302, 3249, 3429, 5580, 5796, 6192, 6732, 7092, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fg(
                        pbuffer, 10902, 3879, 4014, 6012, 6192, 7452, 7542, 7812, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fg(
                        pbuffer, 11352, 4284, 4554, 6372, 6732, 7542, 8082, 8622, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fg(
                        pbuffer, 12252, 4554, 4824, 6732, 7092, 7812, 8622, 9162, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_gg(
                        pbuffer, 13152, 8082, 8622, 9702, 10302, 10902, 11352, 12252, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 13152, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom020SumRecGG_hpp */
