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


#include "SimdElectronRepulsionRecLF.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformLF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lf_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(3108, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 7, 8, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 8, 9, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 9, 10, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 10, 11, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 11, 12, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 12, 13, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 13, 14, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 14, 15, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 15, 16, 45, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 75, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 81, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 87, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 93, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 99, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 105, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 111, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 117, 0, 39, 42, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 123, 0, 48, 51, 81, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 132, 0, 51, 54, 87, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 141, 0, 54, 57, 93, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_4(buffer, 150, 0, 57, 60, 99, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 161, 0, 60, 63, 105, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 170, 0, 63, 66, 111, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 179, 0, 66, 69, 117, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 188, 0, 75, 81, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 200, 0, 81, 87, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_18(buffer, 212, 0, 87, 93, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_5(buffer, 228, 0, 93, 99, 161, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 243, 0, 99, 105, 170, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 256, 0, 105, 111, 179, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 269, 0, 123, 132, 200, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_18(buffer, 284, 0, 132, 141, 212, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_19(buffer, 299, 0, 141, 150, 228, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_14(buffer, 322, 0, 150, 161, 243, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_16(buffer, 342, 0, 161, 170, 256, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 360, 0, 188, 200, 284, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_12(buffer, 375, 0, 200, 212, 299, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_13(buffer, 393, 0, 212, 228, 322, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_14(buffer, 420, 0, 228, 243, 342, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 446, 0, 269, 284, 375, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_4(buffer, 464, 0, 284, 299, 393, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_5(buffer, 482, 0, 299, 322, 420, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 513, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 516, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 519, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 522, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 525, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 528, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 531, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 534, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 537, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 540, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 543, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 546, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 549, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 552, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 555, 3, 15, 42, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 558, 3, 18, 48, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 561, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 564, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 567, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 570, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 573, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 576, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 579, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 582, 3, 42, 72, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 585, 3, 51, 81, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 588, 3, 54, 87, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 591, 3, 57, 93, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 594, 3, 60, 99, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 597, 3, 63, 105, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 600, 3, 66, 111, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 603, 3, 69, 117, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 606, 3, 75, 123, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 609, 3, 81, 132, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 612, 3, 87, 141, ncols, p);

            compute_prim_gp_electron_repulsion_16(buffer, 615, 3, 93, 150, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 618, 3, 99, 161, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 621, 3, 105, 170, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 624, 3, 111, 179, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 627, 3, 132, 200, ncols, p);

            compute_prim_hp_electron_repulsion_18(buffer, 630, 3, 141, 212, ncols, p);

            compute_prim_hp_electron_repulsion_19(buffer, 633, 0, 3, 150, 618, 228, ncols, p);

            compute_prim_hp_electron_repulsion_14(buffer, 649, 3, 161, 243, ncols, p);

            compute_prim_hp_electron_repulsion_11(buffer, 658, 3, 170, 256, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 661, 3, 188, 269, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 664, 3, 200, 284, ncols, p);

            compute_prim_ip_electron_repulsion_16(buffer, 667, 0, 3, 212, 633, 299, ncols, p);

            compute_prim_ip_electron_repulsion_10(buffer, 700, 0, 3, 228, 649, 322, ncols, p);

            compute_prim_ip_electron_repulsion_14(buffer, 727, 3, 243, 342, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 736, 3, 284, 375, ncols, p);

            compute_prim_kp_electron_repulsion_9(buffer, 757, 0, 3, 299, 700, 393, ncols, p);

            compute_prim_kp_electron_repulsion_10(buffer, 806, 0, 3, 322, 727, 420, ncols, p);

            compute_prim_lp_electron_repulsion_2(buffer, 838, 3, 360, 446, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 862, 3, 375, 464, ncols, p);

            compute_prim_lp_electron_repulsion_4(buffer, 886, 0, 3, 393, 806, 482, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 946, 3, 8, 9, 516, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 949, 3, 9, 10, 519, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 952, 3, 10, 11, 522, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 955, 3, 11, 12, 525, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 958, 3, 12, 13, 528, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 961, 3, 13, 14, 531, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 964, 3, 14, 15, 534, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 967, 0, 513, 946, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 970, 0, 516, 949, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 973, 0, 519, 952, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 976, 0, 522, 955, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 979, 0, 525, 958, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 982, 0, 528, 961, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 985, 0, 531, 964, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 988, 3, 537, 48, 51, 564, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 991, 3, 540, 51, 54, 567, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 994, 3, 543, 54, 57, 570, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 997, 3, 546, 57, 60, 573, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1000, 3, 549, 60, 63, 576, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1003, 3, 552, 63, 66, 579, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1006, 3, 555, 66, 69, 582, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1009, 0, 3, 564, 991, 75, 81, 588, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1018, 0, 3, 567, 994, 81, 87, 591, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1027, 0, 3, 570, 997, 87, 93, 594, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1036, 0, 3, 573, 1000, 93, 99, 597, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1045, 0, 3, 576, 1003, 99, 105, 600, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1054, 0, 3, 579, 1006, 105, 111, 603, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1063, 0, 3, 988, 991, 588, 1018, 123, 132, 612, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1078, 0, 3, 991, 994, 591, 1027, 132, 141, 615, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_38(buffer, 1093, 0, 3, 994, 997, 594, 1036, 141, 150, 618, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_35(buffer, 1108, 0, 3, 997, 1000, 597, 1045, 150, 161, 621, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1123, 0, 3, 1000, 1003, 600, 1054, 161, 170, 624, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_22(buffer, 1138, 0, 3, 1009, 1018, 612, 1078, 188, 200, 630, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_38(buffer, 1162, 0, 3, 1018, 1027, 615, 1093, 200, 212, 633, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_39(buffer, 1186, 0, 3, 1027, 1036, 618, 1108, 212, 228, 649, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_40(buffer, 1216, 0, 3, 1036, 1045, 621, 1123, 228, 243, 658, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_25(buffer, 1240, 0, 3, 1063, 1078, 630, 1162, 269, 284, 667, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_26(buffer, 1276, 0, 3, 1078, 1093, 633, 1186, 284, 299, 700, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_27(buffer, 1338, 0, 3, 1093, 1108, 649, 1216, 299, 322, 727, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_10(buffer, 1380, 0, 3, 1138, 1162, 667, 1276, 360, 375, 757, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_11(buffer, 1506, 0, 3, 1162, 1186, 700, 1338, 375, 393, 806, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_1(buffer, 1587, 0, 3, 1240, 1276, 757, 1506, 446, 464, 886, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1764, 3, 967, 558, 561, 988, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1767, 3, 970, 561, 564, 991, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1770, 3, 973, 564, 567, 994, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1773, 3, 976, 567, 570, 997, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1776, 3, 979, 570, 573, 1000, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1779, 3, 982, 573, 576, 1003, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1782, 3, 985, 576, 579, 1006, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1785, 0, 3, 991, 1770, 585, 588, 1018, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1794, 0, 3, 994, 1773, 588, 591, 1027, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1803, 0, 3, 997, 1776, 591, 594, 1036, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1812, 0, 3, 1000, 1779, 594, 597, 1045, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1821, 0, 3, 1003, 1782, 597, 600, 1054, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1830, 0, 3, 1764, 1767, 1009, 1785, 606, 609, 1063, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1845, 0, 3, 1767, 1770, 1018, 1794, 609, 612, 1078, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1860, 0, 3, 1770, 1773, 1027, 1803, 612, 615, 1093, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1875, 0, 3, 1773, 1776, 1036, 1812, 615, 618, 1108, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1890, 0, 3, 1776, 1779, 1045, 1821, 618, 621, 1123, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_23(buffer, 1905, 0, 3, 1785, 1794, 1078, 1860, 627, 630, 1162, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_41(buffer, 1929, 0, 3, 1794, 1803, 1093, 1875, 630, 633, 1186, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_42(buffer, 1953, 0, 3, 1803, 1812, 1108, 1890, 633, 649, 1216, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_23(buffer, 1977, 0, 3, 1830, 1845, 1138, 1905, 661, 664, 1240, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_24(buffer, 2013, 0, 3, 1845, 1860, 1162, 1929, 664, 667, 1276, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_25(buffer, 2049, 0, 3, 1860, 1875, 1186, 1953, 667, 700, 1338, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_10(buffer, 2091, 0, 3, 1905, 1929, 1276, 2049, 736, 757, 1506, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 2208, 0, 3, 1977, 2013, 1380, 2091, 838, 862, 1587, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2658, 2208, 450, ncols);
        }
    }

    simdtrf::transform_lf(values, nvalues, buffer, 2658, nmax);
}

}  // namespace simdt2ceri
