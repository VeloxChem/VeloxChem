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


#include "SimdElectronRepulsionRecIG.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformIG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ig_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ig_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2750, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 13, 14, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 14, 15, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 15, 16, 42, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 72, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 78, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 84, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 90, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 114, 0, 45, 48, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 123, 0, 48, 51, 78, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 132, 0, 51, 54, 84, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 141, 0, 54, 57, 90, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 150, 0, 57, 60, 96, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 159, 0, 60, 63, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 168, 0, 63, 66, 108, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 177, 0, 72, 78, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 189, 0, 78, 84, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 201, 0, 84, 90, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 213, 0, 90, 96, 159, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 225, 0, 96, 102, 168, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 237, 0, 114, 123, 177, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 249, 0, 123, 132, 189, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 261, 0, 132, 141, 201, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 273, 0, 141, 150, 213, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 285, 0, 150, 159, 225, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 297, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 300, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 303, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 306, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 309, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 312, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 315, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 318, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 321, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 324, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 327, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 330, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 333, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 336, 3, 15, 39, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 339, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 342, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 345, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 348, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 357, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 366, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 369, 3, 39, 69, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 372, 3, 51, 78, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 375, 3, 54, 84, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 378, 3, 57, 90, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 387, 3, 60, 96, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 396, 3, 63, 102, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 405, 3, 66, 108, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 408, 3, 78, 132, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 411, 3, 84, 141, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 423, 3, 90, 150, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 435, 3, 96, 159, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 447, 3, 102, 168, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 456, 3, 132, 189, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 471, 3, 141, 201, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 486, 3, 150, 213, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 501, 3, 159, 225, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 516, 3, 189, 261, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 534, 3, 201, 273, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 552, 3, 213, 285, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 570, 3, 9, 10, 300, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 573, 3, 10, 11, 303, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 576, 3, 11, 12, 306, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 579, 3, 12, 13, 309, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 582, 3, 13, 14, 312, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 585, 3, 14, 15, 315, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 588, 0, 300, 573, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 591, 0, 303, 576, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 594, 0, 306, 579, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 597, 0, 309, 582, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 600, 0, 312, 585, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 603, 3, 318, 45, 48, 339, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 606, 3, 321, 48, 51, 342, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 609, 3, 324, 51, 54, 345, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 612, 3, 327, 54, 57, 348, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 615, 3, 330, 57, 60, 357, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 624, 3, 333, 60, 63, 366, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 627, 3, 336, 63, 66, 369, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 630, 0, 3, 342, 609, 72, 78, 375, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 639, 0, 3, 345, 612, 78, 84, 378, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_9(buffer, 648, 0, 3, 348, 615, 84, 90, 387, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 669, 0, 3, 357, 624, 90, 96, 396, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 684, 0, 3, 366, 627, 96, 102, 405, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 693, 0, 3, 603, 606, 372, 630, 114, 123, 408, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 708, 0, 3, 606, 609, 375, 639, 123, 132, 411, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_27(buffer, 723, 0, 3, 609, 612, 378, 648, 132, 141, 423, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_28(buffer, 759, 0, 3, 612, 615, 387, 669, 141, 150, 435, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_29(buffer, 789, 0, 3, 615, 624, 396, 684, 150, 159, 447, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_12(buffer, 810, 0, 3, 630, 639, 411, 723, 177, 189, 471, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_13(buffer, 843, 0, 3, 639, 648, 423, 759, 189, 201, 486, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_14(buffer, 894, 0, 3, 648, 669, 435, 789, 201, 213, 501, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_2(buffer, 936, 0, 3, 693, 708, 456, 810, 237, 249, 516, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_3(buffer, 972, 0, 3, 708, 723, 471, 843, 249, 261, 534, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_4(buffer, 1008, 0, 3, 723, 759, 486, 894, 261, 273, 552, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1077, 3, 297, 300, 573, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1080, 3, 300, 303, 576, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1083, 3, 303, 306, 579, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1086, 3, 306, 309, 582, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1089, 3, 309, 312, 585, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1092, 0, 570, 1077, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1095, 0, 573, 1080, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1098, 0, 576, 1083, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1101, 0, 579, 1086, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1104, 0, 582, 1089, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1107, 3, 588, 339, 342, 609, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1110, 3, 591, 342, 345, 612, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 1113, 3, 594, 345, 348, 615, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 1116, 3, 597, 348, 357, 624, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_23(buffer, 1119, 3, 600, 357, 366, 627, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 1122, 0, 3, 609, 1110, 372, 375, 639, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_28(buffer, 1134, 0, 3, 612, 1113, 375, 378, 648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_29(buffer, 1146, 0, 3, 615, 1116, 378, 387, 669, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_30(buffer, 1164, 0, 3, 624, 1119, 387, 396, 684, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_21(buffer, 1173, 0, 3, 1107, 1110, 639, 1134, 408, 411, 723, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1191, 0, 3, 1110, 1113, 648, 1146, 411, 423, 759, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1246, 0, 3, 1113, 1116, 669, 1164, 423, 435, 789, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_9(buffer, 1279, 0, 3, 1122, 1134, 723, 1191, 456, 471, 843, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_10(buffer, 1417, 0, 3, 1134, 1146, 759, 1246, 471, 486, 894, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_1(buffer, 1499, 0, 3, 1173, 1191, 843, 1417, 516, 534, 1008, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1688, 3, 1092, 603, 606, 1107, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1691, 3, 1095, 606, 609, 1110, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1694, 3, 1098, 609, 612, 1113, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_24(buffer, 1697, 3, 1101, 612, 615, 1116, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_25(buffer, 1700, 3, 1104, 615, 624, 1119, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_20(buffer, 1703, 0, 3, 1110, 1694, 630, 639, 1134, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_21(buffer, 1712, 0, 3, 1113, 1697, 639, 648, 1146, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_22(buffer, 1721, 0, 3, 1116, 1700, 648, 669, 1164, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_11(buffer, 1730, 0, 3, 1688, 1691, 1122, 1703, 693, 708, 1173, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_12(buffer, 1748, 0, 3, 1691, 1694, 1134, 1712, 708, 723, 1191, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_13(buffer, 1766, 0, 3, 1694, 1697, 1146, 1721, 723, 759, 1246, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_5(buffer, 1799, 0, 3, 1703, 1712, 1191, 1766, 810, 843, 1417, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 1910, 0, 3, 1730, 1748, 1279, 1799, 936, 972, 1499, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2330, 1910, 420, ncols);
        }
    }

    simdtrf::transform_ig(values, nvalues, buffer, 2330, nmax);
}

}  // namespace simdt2ceri
