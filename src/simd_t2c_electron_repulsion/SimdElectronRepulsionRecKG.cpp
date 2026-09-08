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


#include "SimdElectronRepulsionRecKG.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformKG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kg_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(3827, nvalues);

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

            compute_prim_hs_electron_repulsion_14(buffer, 195, 0, 78, 84, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 207, 0, 84, 90, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 219, 0, 90, 96, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 231, 0, 96, 102, 159, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 243, 0, 102, 108, 168, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 255, 0, 108, 114, 177, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 267, 0, 114, 120, 186, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 279, 0, 132, 141, 219, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 294, 0, 141, 150, 231, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 309, 0, 150, 159, 243, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 324, 0, 159, 168, 255, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 339, 0, 168, 177, 267, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 354, 0, 195, 207, 279, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 369, 0, 207, 219, 294, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 384, 0, 219, 231, 309, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 399, 0, 231, 243, 324, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 414, 0, 243, 255, 339, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 429, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 432, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 435, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 438, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 441, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 444, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 447, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 450, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 453, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 456, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 459, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 462, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 465, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 468, 3, 15, 45, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 471, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 474, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 477, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 480, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 483, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 486, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 489, 3, 42, 72, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 492, 3, 45, 75, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 495, 3, 54, 90, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 498, 3, 57, 96, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 501, 3, 60, 102, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 504, 3, 63, 108, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 513, 3, 66, 114, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 522, 3, 69, 120, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 525, 3, 72, 126, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 528, 3, 90, 141, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 531, 3, 96, 150, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 534, 3, 102, 159, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 546, 3, 108, 168, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 558, 3, 114, 177, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 567, 3, 120, 186, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 570, 3, 141, 219, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 573, 3, 150, 231, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 588, 3, 159, 243, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 603, 3, 168, 255, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 618, 3, 177, 267, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 627, 3, 219, 294, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 645, 3, 231, 309, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 663, 3, 243, 324, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 681, 3, 255, 339, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 699, 3, 294, 384, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 720, 3, 309, 399, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 741, 3, 324, 414, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 762, 3, 9, 10, 432, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 765, 3, 10, 11, 435, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 768, 3, 11, 12, 438, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 771, 3, 12, 13, 441, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 774, 3, 13, 14, 444, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 777, 3, 14, 15, 447, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 780, 0, 429, 762, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 783, 0, 432, 765, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 786, 0, 435, 768, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 789, 0, 438, 771, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 792, 0, 441, 774, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 795, 0, 444, 777, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 798, 3, 450, 51, 54, 474, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 801, 3, 453, 54, 57, 477, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 804, 3, 456, 57, 60, 480, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 807, 3, 459, 60, 63, 483, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 810, 3, 462, 63, 66, 486, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 813, 3, 465, 66, 69, 489, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 816, 3, 468, 69, 72, 492, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 819, 0, 3, 471, 798, 78, 84, 495, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 828, 0, 3, 474, 801, 84, 90, 498, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 837, 0, 3, 477, 804, 90, 96, 501, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 846, 0, 3, 480, 807, 96, 102, 504, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 855, 0, 3, 483, 810, 102, 108, 513, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 870, 0, 3, 486, 813, 108, 114, 522, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 879, 0, 3, 489, 816, 114, 120, 525, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 888, 0, 3, 798, 801, 498, 837, 132, 141, 531, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 903, 0, 3, 801, 804, 501, 846, 141, 150, 534, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_36(buffer, 918, 0, 3, 804, 807, 504, 855, 150, 159, 546, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_37(buffer, 948, 0, 3, 807, 810, 513, 870, 159, 168, 558, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 969, 0, 3, 810, 813, 522, 879, 168, 177, 567, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 984, 0, 3, 819, 828, 528, 888, 195, 207, 570, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1005, 0, 3, 828, 837, 531, 903, 207, 219, 573, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_30(buffer, 1026, 0, 3, 837, 846, 534, 918, 219, 231, 588, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_31(buffer, 1080, 0, 3, 846, 855, 546, 948, 231, 243, 603, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_32(buffer, 1122, 0, 3, 855, 870, 558, 969, 243, 255, 618, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_12(buffer, 1152, 0, 3, 888, 903, 573, 1026, 279, 294, 645, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_13(buffer, 1194, 0, 3, 903, 918, 588, 1080, 294, 309, 663, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_14(buffer, 1263, 0, 3, 918, 948, 603, 1122, 309, 324, 681, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_2(buffer, 1320, 0, 3, 984, 1005, 627, 1152, 354, 369, 699, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_3(buffer, 1365, 0, 3, 1005, 1026, 645, 1194, 369, 384, 720, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_4(buffer, 1410, 0, 3, 1026, 1080, 663, 1263, 384, 399, 741, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1500, 3, 429, 432, 765, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1503, 3, 432, 435, 768, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1506, 3, 435, 438, 771, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1509, 3, 438, 441, 774, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1512, 3, 441, 444, 777, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1515, 0, 762, 1500, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1518, 0, 765, 1503, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1521, 0, 768, 1506, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1524, 0, 771, 1509, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1527, 0, 774, 1512, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1530, 3, 780, 471, 474, 801, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1533, 3, 783, 474, 477, 804, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1536, 3, 786, 477, 480, 807, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1539, 3, 789, 480, 483, 810, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1542, 3, 792, 483, 486, 813, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1545, 3, 795, 486, 489, 816, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1548, 0, 3, 801, 1533, 495, 498, 837, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 1557, 0, 3, 804, 1536, 498, 501, 846, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_38(buffer, 1569, 0, 3, 807, 1539, 501, 504, 855, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_39(buffer, 1581, 0, 3, 810, 1542, 504, 513, 870, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_40(buffer, 1593, 0, 3, 813, 1545, 513, 522, 879, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_41(buffer, 1602, 0, 3, 1530, 1533, 837, 1557, 528, 531, 903, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_42(buffer, 1623, 0, 3, 1533, 1536, 846, 1569, 531, 534, 918, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_43(buffer, 1644, 0, 3, 1536, 1539, 855, 1581, 534, 546, 948, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_44(buffer, 1671, 0, 3, 1539, 1542, 870, 1593, 546, 558, 969, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_26(buffer, 1689, 0, 3, 1548, 1557, 903, 1623, 570, 573, 1026, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_27(buffer, 1719, 0, 3, 1557, 1569, 918, 1644, 573, 588, 1080, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_28(buffer, 1791, 0, 3, 1569, 1581, 948, 1671, 588, 603, 1122, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_11(buffer, 1836, 0, 3, 1602, 1623, 1026, 1719, 627, 645, 1194, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_12(buffer, 2021, 0, 3, 1623, 1644, 1080, 1791, 645, 663, 1263, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_1(buffer, 2129, 0, 3, 1689, 1719, 1194, 2021, 699, 720, 1410, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2372, 3, 1515, 798, 801, 1533, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2375, 3, 1518, 801, 804, 1536, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2378, 3, 1521, 804, 807, 1539, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2381, 3, 1524, 807, 810, 1542, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2384, 3, 1527, 810, 813, 1545, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 2387, 0, 3, 1530, 2372, 819, 828, 1548, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 2396, 0, 3, 1533, 2375, 828, 837, 1557, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 2405, 0, 3, 1536, 2378, 837, 846, 1569, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_42(buffer, 2414, 0, 3, 1539, 2381, 846, 855, 1581, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_43(buffer, 2423, 0, 3, 1542, 2384, 855, 870, 1593, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_39(buffer, 2432, 0, 3, 2372, 2375, 1557, 2405, 888, 903, 1623, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_40(buffer, 2450, 0, 3, 2375, 2378, 1569, 2414, 903, 918, 1644, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_41(buffer, 2468, 0, 3, 2378, 2381, 1581, 2423, 918, 948, 1671, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_23(buffer, 2486, 0, 3, 2387, 2396, 1602, 2432, 984, 1005, 1689, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_24(buffer, 2519, 0, 3, 2396, 2405, 1623, 2450, 1005, 1026, 1719, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_25(buffer, 2552, 0, 3, 2405, 2414, 1644, 2468, 1026, 1080, 1791, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_10(buffer, 2600, 0, 3, 2432, 2450, 1719, 2552, 1152, 1194, 2021, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 2747, 0, 3, 2486, 2519, 1836, 2600, 1320, 1365, 2129, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3287, 2747, 540, ncols);
        }
    }

    simdtrf::transform_kg(values, nvalues, buffer, 3287, nmax);
}

}  // namespace simdt2ceri
