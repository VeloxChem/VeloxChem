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


#include "SimdElectronRepulsionRecHI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformHI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hi_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(4444, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = nvalues;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto pi = mathconst::pi_value();

            const auto fj = 2.0 * a_norms[i] * b_norms[j] * pi * pi * std::sqrt(pi)
                            / (a_exps[i] * b_exps[j] * std::sqrt(p));

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 8, 9, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 9, 10, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 10, 11, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 11, 12, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 12, 13, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 13, 14, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 14, 15, 45, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 75, 0, 15, 16, 48, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 78, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 84, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 90, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 114, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 120, 0, 39, 42, 72, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 126, 0, 42, 45, 75, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 132, 0, 51, 54, 90, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 141, 0, 54, 57, 96, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 150, 0, 57, 60, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 159, 0, 60, 63, 108, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 168, 0, 63, 66, 114, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 177, 0, 66, 69, 120, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 186, 0, 69, 72, 126, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 195, 0, 78, 84, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 204, 0, 84, 90, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 213, 0, 90, 96, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 222, 0, 96, 102, 159, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 231, 0, 102, 108, 168, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 240, 0, 108, 114, 177, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 249, 0, 114, 120, 186, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 258, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 261, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 264, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 267, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 270, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 273, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 276, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 279, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 282, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 285, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 288, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 291, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 294, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 297, 3, 15, 45, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 300, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 303, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 306, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 315, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 324, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 333, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 342, 3, 42, 72, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 351, 3, 45, 75, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 354, 3, 54, 90, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 357, 3, 57, 96, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 366, 3, 60, 102, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 375, 3, 63, 108, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 384, 3, 66, 114, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 393, 3, 69, 120, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 402, 3, 72, 126, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 411, 3, 90, 141, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 423, 3, 96, 150, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 435, 3, 102, 159, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 447, 3, 108, 168, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 459, 3, 114, 177, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 471, 3, 120, 186, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 483, 3, 141, 213, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 498, 3, 150, 222, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 513, 3, 159, 231, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 528, 3, 168, 240, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 543, 3, 177, 249, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 558, 3, 9, 10, 261, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 561, 3, 10, 11, 264, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 564, 3, 11, 12, 267, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 567, 3, 12, 13, 270, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 570, 3, 13, 14, 273, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 573, 3, 14, 15, 276, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 576, 0, 258, 558, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 579, 0, 261, 561, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 582, 0, 264, 564, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 585, 0, 267, 567, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 588, 0, 270, 570, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 591, 0, 273, 573, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 594, 3, 279, 51, 54, 303, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 597, 3, 282, 54, 57, 306, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 600, 3, 285, 57, 60, 315, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 609, 3, 288, 60, 63, 324, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 618, 3, 291, 63, 66, 333, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 627, 3, 294, 66, 69, 342, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 636, 3, 297, 69, 72, 351, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 639, 0, 3, 300, 594, 78, 84, 354, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 648, 0, 3, 303, 597, 84, 90, 357, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 657, 0, 3, 306, 600, 90, 96, 366, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 672, 0, 3, 315, 609, 96, 102, 375, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 687, 0, 3, 324, 618, 102, 108, 384, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 702, 0, 3, 333, 627, 108, 114, 393, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 717, 0, 3, 342, 636, 114, 120, 402, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 732, 0, 3, 594, 597, 357, 657, 132, 141, 423, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 756, 0, 3, 597, 600, 366, 672, 141, 150, 435, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 780, 0, 3, 600, 609, 375, 687, 150, 159, 447, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 804, 0, 3, 609, 618, 384, 702, 159, 168, 459, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 828, 0, 3, 618, 627, 393, 717, 168, 177, 471, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_2(buffer, 852, 0, 3, 639, 648, 411, 732, 195, 204, 483, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_5(buffer, 879, 0, 3, 648, 657, 423, 756, 204, 213, 498, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 906, 0, 3, 657, 672, 435, 780, 213, 222, 513, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 933, 0, 3, 672, 687, 447, 804, 222, 231, 528, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 960, 0, 3, 687, 702, 459, 828, 231, 240, 543, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 987, 3, 258, 261, 561, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 990, 3, 261, 264, 564, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 993, 3, 264, 267, 567, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 996, 3, 267, 270, 570, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 999, 3, 270, 273, 573, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1002, 0, 558, 987, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1005, 0, 561, 990, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1008, 0, 564, 993, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1011, 0, 567, 996, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1014, 0, 570, 999, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1017, 3, 576, 300, 303, 597, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 1020, 3, 579, 303, 306, 600, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_17(buffer, 1023, 3, 582, 306, 315, 609, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1044, 3, 585, 315, 324, 618, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1062, 3, 588, 324, 333, 627, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 1080, 3, 591, 333, 342, 636, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_20(buffer, 1083, 0, 3, 597, 1020, 354, 357, 657, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_11(buffer, 1092, 0, 3, 600, 1023, 357, 366, 672, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_23(buffer, 1119, 0, 3, 609, 1044, 366, 375, 687, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_24(buffer, 1152, 0, 3, 618, 1062, 375, 384, 702, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_19(buffer, 1179, 0, 3, 627, 1080, 384, 393, 717, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_14(buffer, 1203, 0, 3, 1017, 1020, 657, 1092, 411, 423, 756, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_15(buffer, 1245, 0, 3, 1020, 1023, 672, 1119, 423, 435, 780, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_16(buffer, 1290, 0, 3, 1023, 1044, 687, 1152, 435, 447, 804, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_17(buffer, 1338, 0, 3, 1044, 1062, 702, 1179, 447, 459, 828, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_5(buffer, 1380, 0, 3, 1083, 1092, 756, 1245, 483, 498, 906, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_6(buffer, 1434, 0, 3, 1092, 1119, 780, 1290, 498, 513, 933, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_7(buffer, 1488, 0, 3, 1119, 1152, 804, 1338, 513, 528, 960, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1542, 3, 558, 561, 990, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1545, 3, 561, 564, 993, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1548, 3, 564, 567, 996, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1551, 3, 567, 570, 999, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1554, 0, 987, 1542, ncols, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1557, 0, 990, 1545, ncols, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1560, 0, 993, 1548, ncols, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1563, 0, 996, 1551, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1566, 3, 1002, 594, 597, 1020, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_22(buffer, 1569, 3, 1005, 597, 600, 1023, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_23(buffer, 1572, 0, 3, 1008, 1560, 600, 609, 1044, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1602, 3, 1011, 609, 618, 1062, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_21(buffer, 1629, 3, 1014, 618, 627, 1080, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_10(buffer, 1632, 0, 3, 1017, 1566, 639, 648, 1083, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 1641, 0, 3, 1020, 1569, 648, 657, 1092, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_17(buffer, 1650, 0, 3, 1023, 1572, 657, 672, 1119, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_18(buffer, 1707, 0, 3, 1044, 1602, 672, 687, 1152, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_19(buffer, 1755, 0, 3, 1062, 1629, 687, 702, 1179, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_8(buffer, 1788, 0, 3, 1566, 1569, 1092, 1650, 732, 756, 1245, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_9(buffer, 1848, 0, 3, 1569, 1572, 1119, 1707, 756, 780, 1290, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_10(buffer, 1937, 0, 3, 1572, 1602, 1152, 1755, 780, 804, 1338, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_2(buffer, 2003, 0, 3, 1632, 1641, 1203, 1788, 852, 879, 1380, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_3(buffer, 2084, 0, 3, 1641, 1650, 1245, 1848, 879, 906, 1434, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_4(buffer, 2165, 0, 3, 1650, 1707, 1290, 1937, 906, 933, 1488, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_8(buffer, 2270, 3, 987, 990, 1545, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_8(buffer, 2273, 3, 990, 993, 1548, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_8(buffer, 2276, 3, 993, 996, 1551, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_8(buffer, 2279, 0, 1542, 2270, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 2282, 0, 1545, 2273, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 2285, 0, 1548, 2276, ncols, p);

            compute_prim_dh_electron_repulsion_5(buffer, 2288, 3, 1554, 1017, 1020, 1569, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_13(buffer, 2291, 3, 1557, 1020, 1023, 1572, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_14(buffer, 2294, 3, 1560, 1023, 1044, 1602, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_15(buffer, 2321, 3, 1563, 1044, 1062, 1629, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_8(buffer, 2324, 0, 3, 1569, 2291, 1083, 1092, 1650, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_9(buffer, 2333, 0, 3, 1572, 2294, 1092, 1119, 1707, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_10(buffer, 2413, 0, 3, 1602, 2321, 1119, 1152, 1755, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_3(buffer, 2458, 0, 3, 2288, 2291, 1650, 2333, 1203, 1245, 1848, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_4(buffer, 2657, 0, 3, 2291, 2294, 1707, 2413, 1245, 1290, 1937, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_1(buffer, 2776, 0, 3, 2324, 2333, 1848, 2657, 1380, 1434, 2165, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_2(buffer, 3052, 3, 2279, 1566, 1569, 2291, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_5(buffer, 3055, 3, 2282, 1569, 1572, 2294, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_6(buffer, 3058, 3, 2285, 1572, 1602, 2321, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_2(buffer, 3061, 0, 3, 2288, 3052, 1632, 1641, 2324, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_3(buffer, 3070, 0, 3, 2291, 3055, 1641, 1650, 2333, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_4(buffer, 3079, 0, 3, 2294, 3058, 1650, 1707, 2413, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_1(buffer, 3124, 0, 3, 3052, 3055, 2333, 3079, 1788, 1848, 2657, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 3268, 0, 3, 3061, 3070, 2458, 3124, 2003, 2084, 2776, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3856, 3268, 588, ncols);
        }
    }

    simdtrf::transform_hi(values, nvalues, buffer, 3856, nmax);
}

}  // namespace simdt2ceri
