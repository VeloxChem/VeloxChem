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

#ifndef ElectricDipoleMomentumGeom110RecGG_hpp
#define ElectricDipoleMomentumGeom110RecGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "ElectricDipoleMomentumPrimRecDD.hpp"
#include "ElectricDipoleMomentumPrimRecDF.hpp"
#include "ElectricDipoleMomentumPrimRecDG.hpp"
#include "ElectricDipoleMomentumPrimRecDP.hpp"
#include "ElectricDipoleMomentumPrimRecFD.hpp"
#include "ElectricDipoleMomentumPrimRecFF.hpp"
#include "ElectricDipoleMomentumPrimRecFG.hpp"
#include "ElectricDipoleMomentumPrimRecGF.hpp"
#include "ElectricDipoleMomentumPrimRecGG.hpp"
#include "ElectricDipoleMomentumPrimRecHG.hpp"
#include "ElectricDipoleMomentumPrimRecPD.hpp"
#include "ElectricDipoleMomentumPrimRecPF.hpp"
#include "ElectricDipoleMomentumPrimRecPG.hpp"
#include "ElectricDipoleMomentumPrimRecPP.hpp"
#include "ElectricDipoleMomentumPrimRecPS.hpp"
#include "ElectricDipoleMomentumPrimRecSD.hpp"
#include "ElectricDipoleMomentumPrimRecSF.hpp"
#include "ElectricDipoleMomentumPrimRecSG.hpp"
#include "ElectricDipoleMomentumPrimRecSP.hpp"
#include "ElectricDipoleMomentumPrimRecSS.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"
#include "GtoBlock.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "OverlapPrimRecDG.hpp"
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFG.hpp"
#include "OverlapPrimRecGG.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "OverlapPrimRecSG.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes (d^(1)/dA^(1)G|r|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electric_dipole_momentum_geom_10_gg(T&                               distributor,
                                         const CGtoBlock&                 bra_gto_block,
                                         const CGtoBlock&                 ket_gto_block,
                                         const std::pair<size_t, size_t>& bra_indices,
                                         const std::pair<size_t, size_t>& ket_indices,
                                         const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto r_c = distributor.coordinates()[0];

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

    CSimdArray<double> pbuffer(6855, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(2025, 1);

    CSimdArray<double> sbuffer(729, 1);

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

                t2cfunc::comp_distances_pc(factors, 17, 8, r_c);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                diprec::comp_prim_electric_dipole_momentum_ss(pbuffer, 1, 0, factors, 17);

                ovlrec::comp_prim_overlap_sp(pbuffer, 4, 0, factors, 14);

                diprec::comp_prim_electric_dipole_momentum_sp(pbuffer, 7, 0, 1, factors, 14, a_exp);

                ovlrec::comp_prim_overlap_sd(pbuffer, 16, 0, 4, factors, 14, a_exp);

                diprec::comp_prim_electric_dipole_momentum_sd(pbuffer, 22, 1, 4, 7, factors, 14, a_exp);

                ovlrec::comp_prim_overlap_sf(pbuffer, 40, 4, 16, factors, 14, a_exp);

                diprec::comp_prim_electric_dipole_momentum_sf(pbuffer, 50, 7, 16, 22, factors, 14, a_exp);

                ovlrec::comp_prim_overlap_sg(pbuffer, 80, 16, 40, factors, 14, a_exp);

                diprec::comp_prim_electric_dipole_momentum_sg(pbuffer, 95, 22, 40, 50, factors, 14, a_exp);

                diprec::comp_prim_electric_dipole_momentum_ps(pbuffer, 140, 0, 1, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 149, 0, 4, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pp(pbuffer, 158, 1, 4, 7, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 185, 4, 16, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pd(pbuffer, 203, 7, 16, 22, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 257, 16, 40, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pf(pbuffer, 287, 22, 40, 50, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_pg(pbuffer, 377, 40, 80, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_pg(pbuffer, 422, 50, 80, 95, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_dp(pbuffer, 557, 7, 140, 149, 158, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 611, 16, 149, 185, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_dd(pbuffer, 647, 22, 158, 185, 203, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_df(pbuffer, 755, 40, 185, 257, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_df(pbuffer, 815, 50, 203, 257, 287, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_dg(pbuffer, 995, 80, 257, 377, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_dg(pbuffer, 1085, 95, 287, 377, 422, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_fd(pbuffer, 1355, 203, 557, 611, 647, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ff(pbuffer, 1535, 257, 611, 755, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_ff(pbuffer, 1635, 287, 647, 755, 815, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_fg(pbuffer, 1935, 377, 755, 995, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_fg(pbuffer, 2085, 422, 815, 995, 1085, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_gf(pbuffer, 2535, 815, 1355, 1535, 1635, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_gg(pbuffer, 2985, 995, 1535, 1935, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_gg(pbuffer, 3210, 1085, 1635, 1935, 2085, factors, 11, a_exp);

                diprec::comp_prim_electric_dipole_momentum_hg(pbuffer, 3885, 2085, 2535, 2985, 3210, factors, 11, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 4830, 2085, 3885, 3, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 4830, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace diprec

#endif /* ElectricDipoleMomentumGeom110RecGG_hpp */
