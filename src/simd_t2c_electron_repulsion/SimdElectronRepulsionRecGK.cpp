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


#include "SimdElectronRepulsionRecGK.hpp"

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
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(4209, nvalues);

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

            compute_prim_fs_electron_repulsion_1(buffer, 75, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 78, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 84, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 90, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 114, 0, 39, 42, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 120, 0, 48, 51, 78, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 126, 0, 51, 54, 84, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 132, 0, 54, 57, 90, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 138, 0, 57, 60, 96, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 144, 0, 60, 63, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 150, 0, 63, 66, 108, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 156, 0, 66, 69, 114, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 162, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 165, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 168, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 171, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 174, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 177, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 180, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 183, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 186, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 189, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 192, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 195, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 198, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 201, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 204, 3, 15, 42, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 207, 3, 18, 48, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 210, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 213, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 222, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 231, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 240, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 249, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 258, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 267, 3, 42, 72, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 276, 3, 51, 78, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 285, 3, 54, 84, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 294, 3, 57, 90, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 303, 3, 60, 96, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 312, 3, 63, 102, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 321, 3, 66, 108, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 330, 3, 69, 114, ncols, p);

            compute_prim_gp_electron_repulsion_2(buffer, 339, 3, 75, 120, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 351, 3, 78, 126, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 363, 3, 84, 132, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 375, 3, 90, 138, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 387, 3, 96, 144, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 399, 3, 102, 150, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 411, 3, 108, 156, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 423, 3, 8, 9, 165, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 426, 3, 9, 10, 168, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 429, 3, 10, 11, 171, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 432, 3, 11, 12, 174, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 435, 3, 12, 13, 177, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 438, 3, 13, 14, 180, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 441, 3, 14, 15, 183, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 444, 0, 162, 423, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 447, 0, 165, 426, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 450, 0, 168, 429, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 453, 0, 171, 432, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 456, 0, 174, 435, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 459, 0, 177, 438, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 462, 0, 180, 441, ncols, p);

            compute_prim_dd_electron_repulsion_11(buffer, 465, 3, 186, 48, 51, 213, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 468, 3, 189, 51, 54, 222, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 477, 3, 192, 54, 57, 231, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 486, 3, 195, 57, 60, 240, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 495, 3, 198, 60, 63, 249, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 504, 3, 201, 63, 66, 258, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 513, 3, 204, 66, 69, 267, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_10(buffer, 522, 0, 3, 213, 468, 75, 78, 285, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 534, 0, 3, 222, 477, 78, 84, 294, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 549, 0, 3, 231, 486, 84, 90, 303, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 564, 0, 3, 240, 495, 90, 96, 312, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 579, 0, 3, 249, 504, 96, 102, 321, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 594, 0, 3, 258, 513, 102, 108, 330, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_5(buffer, 609, 0, 3, 465, 468, 285, 534, 120, 126, 363, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 627, 0, 3, 468, 477, 294, 549, 126, 132, 375, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 645, 0, 3, 477, 486, 303, 564, 132, 138, 387, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 663, 0, 3, 486, 495, 312, 579, 138, 144, 399, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 681, 0, 3, 495, 504, 321, 594, 144, 150, 411, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 699, 3, 162, 165, 426, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 702, 3, 165, 168, 429, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 705, 3, 168, 171, 432, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 711, 3, 171, 174, 435, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 717, 3, 174, 177, 438, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 723, 3, 177, 180, 441, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 729, 0, 426, 702, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 732, 0, 429, 705, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 735, 0, 432, 711, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 738, 0, 435, 717, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 741, 0, 438, 723, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 744, 3, 444, 207, 210, 465, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 747, 3, 447, 210, 213, 468, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 750, 3, 450, 213, 222, 477, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 768, 3, 453, 222, 231, 486, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 786, 3, 456, 231, 240, 495, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 804, 3, 459, 240, 249, 504, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 822, 3, 462, 249, 258, 513, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_13(buffer, 840, 0, 3, 468, 750, 276, 285, 534, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 864, 0, 3, 477, 768, 285, 294, 549, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 891, 0, 3, 486, 786, 294, 303, 564, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 918, 0, 3, 495, 804, 303, 312, 579, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_13(buffer, 945, 0, 3, 504, 822, 312, 321, 594, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_2(buffer, 969, 0, 3, 744, 747, 522, 840, 339, 351, 609, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_5(buffer, 1005, 0, 3, 747, 750, 534, 864, 351, 363, 627, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_8(buffer, 1041, 0, 3, 750, 768, 549, 891, 363, 375, 645, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_8(buffer, 1077, 0, 3, 768, 786, 564, 918, 375, 387, 663, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_9(buffer, 1113, 0, 3, 786, 804, 579, 945, 387, 399, 681, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1149, 3, 423, 426, 702, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_8(buffer, 1152, 3, 426, 429, 705, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_4(buffer, 1155, 3, 429, 432, 711, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1167, 3, 432, 435, 717, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1176, 3, 435, 438, 723, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1185, 0, 699, 1149, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1188, 0, 702, 1152, ncols, p);

            compute_prim_pg_electron_repulsion_9(buffer, 1191, 0, 705, 1155, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1194, 0, 711, 1167, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1197, 0, 717, 1176, ncols, p);

            compute_prim_dg_electron_repulsion_14(buffer, 1200, 3, 729, 465, 468, 750, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1203, 3, 732, 468, 477, 768, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_17(buffer, 1230, 3, 735, 477, 486, 786, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1260, 3, 738, 486, 495, 804, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1287, 3, 741, 495, 504, 822, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_10(buffer, 1314, 0, 3, 750, 1203, 522, 534, 864, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_11(buffer, 1350, 0, 3, 768, 1230, 534, 549, 891, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_12(buffer, 1389, 0, 3, 786, 1260, 549, 564, 918, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_13(buffer, 1431, 0, 3, 804, 1287, 564, 579, 945, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_5(buffer, 1464, 0, 3, 1200, 1203, 864, 1350, 609, 627, 1041, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_6(buffer, 1518, 0, 3, 1203, 1230, 891, 1389, 627, 645, 1077, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_7(buffer, 1572, 0, 3, 1230, 1260, 918, 1431, 645, 663, 1113, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 1626, 3, 699, 702, 1152, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_15(buffer, 1629, 3, 702, 705, 1155, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_16(buffer, 1632, 3, 705, 711, 1167, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_17(buffer, 1644, 3, 711, 717, 1176, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 1656, 0, 1152, 1629, ncols, p);

            compute_prim_ph_electron_repulsion_18(buffer, 1659, 0, 3, 1155, 1632, 1194, ncols, p);

            compute_prim_ph_electron_repulsion_19(buffer, 1674, 0, 1167, 1644, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 1677, 3, 1185, 744, 747, 1200, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_15(buffer, 1680, 3, 1188, 747, 750, 1203, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_16(buffer, 1683, 0, 3, 1191, 1659, 750, 768, 1230, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_17(buffer, 1731, 0, 3, 1194, 1674, 768, 786, 1260, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 1776, 3, 1197, 786, 804, 1287, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_7(buffer, 1815, 0, 3, 1203, 1683, 840, 864, 1350, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_8(buffer, 1860, 0, 3, 1230, 1731, 864, 891, 1389, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_9(buffer, 1939, 0, 3, 1260, 1776, 891, 918, 1431, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_2(buffer, 1990, 0, 3, 1677, 1680, 1314, 1815, 969, 1005, 1464, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_3(buffer, 2068, 0, 3, 1680, 1683, 1350, 1860, 1005, 1041, 1518, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_4(buffer, 2146, 0, 3, 1683, 1731, 1389, 1939, 1041, 1077, 1572, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 2242, 3, 1149, 1152, 1629, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_12(buffer, 2245, 3, 1152, 1155, 1632, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_13(buffer, 2248, 3, 1155, 1167, 1644, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 2260, 0, 1626, 2242, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 2263, 0, 1629, 2245, ncols, p);

            compute_prim_pi_electron_repulsion_15(buffer, 2266, 0, 1632, 2248, ncols, p);

            compute_prim_di_electron_repulsion_12(buffer, 2278, 3, 1656, 1200, 1203, 1683, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_13(buffer, 2281, 0, 3, 1659, 2266, 1203, 1230, 1731, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_14(buffer, 2356, 3, 1674, 1230, 1260, 1776, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_5(buffer, 2410, 0, 3, 1683, 2281, 1314, 1350, 1860, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_6(buffer, 2586, 0, 3, 1731, 2356, 1350, 1389, 1939, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_1(buffer, 2688, 0, 3, 2278, 2281, 1860, 2586, 1464, 1518, 2146, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_5(buffer, 2940, 3, 2260, 1677, 1680, 2278, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_6(buffer, 2943, 3, 2263, 1680, 1683, 2281, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_7(buffer, 2946, 3, 2266, 1683, 1731, 2356, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_2(buffer, 3000, 0, 3, 2281, 2946, 1815, 1860, 2586, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 3129, 0, 3, 2940, 2943, 2410, 3000, 1990, 2068, 2688, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3669, 3129, 540, ncols);
        }
    }

    simdtrf::transform_gk(values, nvalues, buffer, 3669, nmax);
}

}  // namespace simdt2ceri
