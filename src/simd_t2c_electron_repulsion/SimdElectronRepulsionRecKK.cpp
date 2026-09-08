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


#include "SimdElectronRepulsionRecKK.hpp"

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
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
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
#include "SimdTransformKK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(13240, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, ncols, fj, mu);

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

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 20, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 7, 8, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 8, 9, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 9, 10, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 10, 11, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 11, 12, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 75, 0, 12, 13, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 78, 0, 13, 14, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 81, 0, 14, 15, 45, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 84, 0, 15, 16, 48, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 87, 0, 16, 17, 51, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 90, 0, 17, 18, 54, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 93, 0, 18, 19, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 21, 24, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 24, 27, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 27, 30, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 114, 0, 30, 33, 72, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 120, 0, 33, 36, 75, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 126, 0, 36, 39, 78, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 132, 0, 39, 42, 81, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 138, 0, 42, 45, 84, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 144, 0, 45, 48, 87, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 150, 0, 48, 51, 90, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 156, 0, 51, 54, 93, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 162, 0, 60, 63, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 171, 0, 63, 66, 108, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 180, 0, 66, 69, 114, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 189, 0, 69, 72, 120, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 198, 0, 72, 75, 126, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 207, 0, 75, 78, 132, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 216, 0, 78, 81, 138, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 225, 0, 81, 84, 144, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 234, 0, 84, 87, 150, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 243, 0, 87, 90, 156, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 252, 0, 96, 102, 171, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 264, 0, 102, 108, 180, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 276, 0, 108, 114, 189, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 288, 0, 114, 120, 198, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 300, 0, 120, 126, 207, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 312, 0, 126, 132, 216, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 324, 0, 132, 138, 225, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 336, 0, 138, 144, 234, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 348, 0, 144, 150, 243, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 360, 0, 162, 171, 264, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 372, 0, 171, 180, 276, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 387, 0, 180, 189, 288, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 402, 0, 189, 198, 300, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 417, 0, 198, 207, 312, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 432, 0, 207, 216, 324, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 447, 0, 216, 225, 336, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 462, 0, 225, 234, 348, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 477, 0, 252, 264, 372, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 492, 0, 264, 276, 387, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 507, 0, 276, 288, 402, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 522, 0, 288, 300, 417, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 537, 0, 300, 312, 432, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 552, 0, 312, 324, 447, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 567, 0, 324, 336, 462, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 582, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 585, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 588, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 591, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 594, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 597, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 600, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 603, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 606, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 609, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 612, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 615, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 618, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 621, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 624, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 627, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 630, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 633, 3, 16, 48, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 636, 3, 17, 51, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 639, 3, 18, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 642, 3, 24, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 645, 3, 27, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 648, 3, 30, 69, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 651, 3, 33, 72, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 654, 3, 36, 75, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 657, 3, 39, 78, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 660, 3, 42, 81, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 663, 3, 45, 84, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 666, 3, 48, 87, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 669, 3, 51, 90, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 672, 3, 54, 93, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 675, 3, 60, 96, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 678, 3, 63, 102, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 681, 3, 66, 108, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 684, 3, 69, 114, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 687, 3, 72, 120, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 696, 3, 75, 126, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 705, 3, 78, 132, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 714, 3, 81, 138, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 723, 3, 84, 144, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 732, 3, 87, 150, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 735, 3, 90, 156, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 738, 3, 102, 171, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 741, 3, 108, 180, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 744, 3, 114, 189, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 756, 3, 120, 198, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 768, 3, 126, 207, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 780, 3, 132, 216, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 792, 3, 138, 225, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 804, 3, 144, 234, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 813, 3, 150, 243, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 816, 3, 162, 252, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 819, 3, 171, 264, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 822, 3, 180, 276, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 837, 3, 189, 288, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 852, 3, 198, 300, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 867, 3, 207, 312, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 882, 3, 216, 324, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 897, 3, 225, 336, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 912, 3, 234, 348, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 921, 3, 264, 372, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 939, 3, 276, 387, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 957, 3, 288, 402, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 975, 3, 300, 417, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 993, 3, 312, 432, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1011, 3, 324, 447, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1029, 3, 336, 462, ncols, p);

            compute_prim_kp_electron_repulsion_2(buffer, 1047, 3, 360, 477, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1068, 3, 372, 492, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1089, 3, 387, 507, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1110, 3, 402, 522, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1131, 3, 417, 537, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1152, 3, 432, 552, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1173, 3, 447, 567, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1194, 3, 9, 10, 585, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1197, 3, 10, 11, 588, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1200, 3, 11, 12, 591, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1203, 3, 12, 13, 594, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1206, 3, 13, 14, 597, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1209, 3, 14, 15, 600, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1212, 3, 15, 16, 603, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1215, 3, 16, 17, 606, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1218, 3, 17, 18, 609, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1221, 0, 582, 1194, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1224, 0, 585, 1197, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1227, 0, 588, 1200, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1230, 0, 591, 1203, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1233, 0, 594, 1206, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1236, 0, 597, 1209, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1239, 0, 600, 1212, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1242, 0, 603, 1215, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1245, 0, 606, 1218, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1248, 3, 612, 60, 63, 645, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1251, 3, 615, 63, 66, 648, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1254, 3, 618, 66, 69, 651, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1257, 3, 621, 69, 72, 654, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1260, 3, 624, 72, 75, 657, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1263, 3, 627, 75, 78, 660, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1266, 3, 630, 78, 81, 663, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1269, 3, 633, 81, 84, 666, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1272, 3, 636, 84, 87, 669, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1275, 3, 639, 87, 90, 672, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1278, 0, 3, 645, 1251, 96, 102, 681, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1287, 0, 3, 648, 1254, 102, 108, 684, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 1296, 0, 3, 651, 1257, 108, 114, 687, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1305, 0, 3, 654, 1260, 114, 120, 696, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1320, 0, 3, 657, 1263, 120, 126, 705, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1335, 0, 3, 660, 1266, 126, 132, 714, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1350, 0, 3, 663, 1269, 132, 138, 723, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1365, 0, 3, 666, 1272, 138, 144, 732, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1374, 0, 3, 669, 1275, 144, 150, 735, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1383, 0, 3, 1248, 1251, 681, 1287, 162, 171, 741, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1398, 0, 3, 1251, 1254, 684, 1296, 171, 180, 744, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1413, 0, 3, 1254, 1257, 687, 1305, 180, 189, 756, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1437, 0, 3, 1257, 1260, 696, 1320, 189, 198, 768, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1461, 0, 3, 1260, 1263, 705, 1335, 198, 207, 780, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1485, 0, 3, 1263, 1266, 714, 1350, 207, 216, 792, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_37(buffer, 1509, 0, 3, 1266, 1269, 723, 1365, 216, 225, 804, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1530, 0, 3, 1269, 1272, 732, 1374, 225, 234, 813, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1545, 0, 3, 1278, 1287, 741, 1398, 252, 264, 822, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1566, 0, 3, 1287, 1296, 744, 1413, 264, 276, 837, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1599, 0, 3, 1296, 1305, 756, 1437, 276, 288, 852, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1632, 0, 3, 1305, 1320, 768, 1461, 288, 300, 867, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1665, 0, 3, 1320, 1335, 780, 1485, 300, 312, 882, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1698, 0, 3, 1335, 1350, 792, 1509, 312, 324, 897, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_34(buffer, 1731, 0, 3, 1350, 1365, 804, 1530, 324, 336, 912, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_15(buffer, 1758, 0, 3, 1383, 1398, 822, 1566, 360, 372, 939, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_20(buffer, 1797, 0, 3, 1398, 1413, 837, 1599, 372, 387, 957, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 1839, 0, 3, 1413, 1437, 852, 1632, 387, 402, 975, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 1881, 0, 3, 1437, 1461, 867, 1665, 402, 417, 993, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 1923, 0, 3, 1461, 1485, 882, 1698, 417, 432, 1011, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_18(buffer, 1965, 0, 3, 1485, 1509, 897, 1731, 432, 447, 1029, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_5(buffer, 2007, 0, 3, 1545, 1566, 939, 1797, 477, 492, 1089, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2052, 0, 3, 1566, 1599, 957, 1839, 492, 507, 1110, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2097, 0, 3, 1599, 1632, 975, 1881, 507, 522, 1131, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2142, 0, 3, 1632, 1665, 993, 1923, 522, 537, 1152, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2187, 0, 3, 1665, 1698, 1011, 1965, 537, 552, 1173, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2232, 3, 582, 585, 1197, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2235, 3, 585, 588, 1200, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2238, 3, 588, 591, 1203, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2241, 3, 591, 594, 1206, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2244, 3, 594, 597, 1209, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2247, 3, 597, 600, 1212, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2250, 3, 600, 603, 1215, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2253, 3, 603, 606, 1218, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2256, 0, 1194, 2232, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2259, 0, 1197, 2235, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2262, 0, 1200, 2238, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2265, 0, 1203, 2241, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2268, 0, 1206, 2244, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2271, 0, 1209, 2247, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2274, 0, 1212, 2250, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2277, 0, 1215, 2253, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2280, 3, 1221, 642, 645, 1251, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2283, 3, 1224, 645, 648, 1254, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2286, 3, 1227, 648, 651, 1257, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2292, 3, 1230, 651, 654, 1260, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2298, 3, 1233, 654, 657, 1263, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2304, 3, 1236, 657, 660, 1266, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2310, 3, 1239, 660, 663, 1269, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2316, 3, 1242, 663, 666, 1272, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2322, 3, 1245, 666, 669, 1275, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2325, 0, 3, 1248, 2280, 675, 678, 1278, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2334, 0, 3, 1251, 2283, 678, 681, 1287, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2343, 0, 3, 1254, 2286, 681, 684, 1296, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_42(buffer, 2358, 0, 3, 1257, 2292, 684, 687, 1305, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2373, 0, 3, 1260, 2298, 687, 696, 1320, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2394, 0, 3, 1263, 2304, 696, 705, 1335, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2415, 0, 3, 1266, 2310, 705, 714, 1350, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_44(buffer, 2436, 0, 3, 1269, 2316, 714, 723, 1365, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_40(buffer, 2451, 0, 3, 1272, 2322, 723, 732, 1374, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_45(buffer, 2460, 0, 3, 2280, 2283, 1287, 2343, 738, 741, 1398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_46(buffer, 2478, 0, 3, 2283, 2286, 1296, 2358, 741, 744, 1413, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2499, 0, 3, 2286, 2292, 1305, 2373, 744, 756, 1437, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2538, 0, 3, 2292, 2298, 1320, 2394, 756, 768, 1461, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2577, 0, 3, 2298, 2304, 1335, 2415, 768, 780, 1485, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_48(buffer, 2616, 0, 3, 2304, 2310, 1350, 2436, 780, 792, 1509, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_49(buffer, 2652, 0, 3, 2310, 2316, 1365, 2451, 792, 804, 1530, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_29(buffer, 2670, 0, 3, 2325, 2334, 1383, 2460, 816, 819, 1545, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_30(buffer, 2697, 0, 3, 2334, 2343, 1398, 2478, 819, 822, 1566, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_38(buffer, 2724, 0, 3, 2343, 2358, 1413, 2499, 822, 837, 1599, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_39(buffer, 2787, 0, 3, 2358, 2373, 1437, 2538, 837, 852, 1632, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_40(buffer, 2850, 0, 3, 2373, 2394, 1461, 2577, 852, 867, 1665, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_36(buffer, 2913, 0, 3, 2394, 2415, 1485, 2616, 867, 882, 1698, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_37(buffer, 2976, 0, 3, 2415, 2436, 1509, 2652, 882, 897, 1731, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_20(buffer, 3018, 0, 3, 2460, 2478, 1566, 2724, 921, 939, 1797, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_21(buffer, 3096, 0, 3, 2478, 2499, 1599, 2787, 939, 957, 1839, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_22(buffer, 3177, 0, 3, 2499, 2538, 1632, 2850, 957, 975, 1881, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_22(buffer, 3258, 0, 3, 2538, 2577, 1665, 2913, 975, 993, 1923, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_19(buffer, 3339, 0, 3, 2577, 2616, 1698, 2976, 993, 1011, 1965, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_2(buffer, 3417, 0, 3, 2670, 2697, 1758, 3018, 1047, 1068, 2007, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_5(buffer, 3507, 0, 3, 2697, 2724, 1797, 3096, 1068, 1089, 2052, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_8(buffer, 3597, 0, 3, 2724, 2787, 1839, 3177, 1089, 1110, 2097, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_8(buffer, 3687, 0, 3, 2787, 2850, 1881, 3258, 1110, 1131, 2142, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_9(buffer, 3777, 0, 3, 2850, 2913, 1923, 3339, 1131, 1152, 2187, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3867, 3, 1194, 1197, 2235, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3870, 3, 1197, 1200, 2238, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3873, 3, 1200, 1203, 2241, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3876, 3, 1203, 1206, 2244, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3879, 3, 1206, 1209, 2247, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3882, 3, 1209, 1212, 2250, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3885, 3, 1212, 1215, 2253, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3888, 0, 2232, 3867, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3891, 0, 2235, 3870, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3894, 0, 2238, 3873, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3897, 0, 2241, 3876, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3900, 0, 2244, 3879, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3903, 0, 2247, 3882, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3906, 0, 2250, 3885, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3909, 3, 2256, 1248, 1251, 2283, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 3912, 3, 2259, 1251, 1254, 2286, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3915, 3, 2262, 1254, 1257, 2292, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3924, 3, 2265, 1257, 1260, 2298, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3933, 3, 2268, 1260, 1263, 2304, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3942, 3, 2271, 1263, 1266, 2310, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3951, 3, 2274, 1266, 1269, 2316, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3960, 3, 2277, 1269, 1272, 2322, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_44(buffer, 3963, 0, 3, 2283, 3912, 1278, 1287, 2343, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 3972, 0, 3, 2286, 3915, 1287, 1296, 2358, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_50(buffer, 3993, 0, 3, 2292, 3924, 1296, 1305, 2373, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 4014, 0, 3, 2298, 3933, 1305, 1320, 2394, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 4041, 0, 3, 2304, 3942, 1320, 1335, 2415, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_52(buffer, 4068, 0, 3, 2310, 3951, 1335, 1350, 2436, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_48(buffer, 4089, 0, 3, 2316, 3960, 1350, 1365, 2451, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_51(buffer, 4098, 0, 3, 3909, 3912, 2343, 3972, 1383, 1398, 2478, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_52(buffer, 4125, 0, 3, 3912, 3915, 2358, 3993, 1398, 1413, 2499, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_53(buffer, 4152, 0, 3, 3915, 3924, 2373, 4014, 1413, 1437, 2538, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_54(buffer, 4209, 0, 3, 3924, 3933, 2394, 4041, 1437, 1461, 2577, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_49(buffer, 4263, 0, 3, 3933, 3942, 2415, 4068, 1461, 1485, 2616, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_50(buffer, 4314, 0, 3, 3942, 3951, 2436, 4089, 1485, 1509, 2652, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_30(buffer, 4335, 0, 3, 3963, 3972, 2478, 4125, 1545, 1566, 2724, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_34(buffer, 4369, 0, 3, 3972, 3993, 2499, 4152, 1566, 1599, 2787, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_35(buffer, 4463, 0, 3, 3993, 4014, 2538, 4209, 1599, 1632, 2850, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_36(buffer, 4563, 0, 3, 4014, 4041, 2577, 4263, 1632, 1665, 2913, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_37(buffer, 4657, 0, 3, 4041, 4068, 2616, 4314, 1665, 1698, 2976, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_16(buffer, 4715, 0, 3, 4098, 4125, 2724, 4369, 1758, 1797, 3096, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_17(buffer, 4835, 0, 3, 4125, 4152, 2787, 4463, 1797, 1839, 3177, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_18(buffer, 4958, 0, 3, 4152, 4209, 2850, 4563, 1839, 1881, 3258, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_19(buffer, 5084, 0, 3, 4209, 4263, 2913, 4657, 1881, 1923, 3339, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_5(buffer, 5201, 0, 3, 4335, 4369, 3096, 4835, 2007, 2052, 3597, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_6(buffer, 5339, 0, 3, 4369, 4463, 3177, 4958, 2052, 2097, 3687, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_7(buffer, 5477, 0, 3, 4463, 4563, 3258, 5084, 2097, 2142, 3777, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5615, 3, 2232, 2235, 3870, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5618, 3, 2235, 2238, 3873, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5621, 3, 2238, 2241, 3876, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5624, 3, 2241, 2244, 3879, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5627, 3, 2244, 2247, 3882, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 5630, 3, 2247, 2250, 3885, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5633, 0, 3867, 5615, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5636, 0, 3870, 5618, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5639, 0, 3873, 5621, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5642, 0, 3876, 5624, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5645, 0, 3879, 5627, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 5648, 0, 3882, 5630, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 5651, 3, 3888, 2280, 2283, 3912, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_32(buffer, 5654, 3, 3891, 2283, 2286, 3915, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 5657, 3, 3894, 2286, 2292, 3924, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 5669, 3, 3897, 2292, 2298, 3933, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 5681, 3, 3900, 2298, 2304, 3942, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 5693, 3, 3903, 2304, 2310, 3951, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 5705, 3, 3906, 2310, 2316, 3960, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_14(buffer, 5708, 0, 3, 3909, 5651, 2325, 2334, 3963, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_46(buffer, 5717, 0, 3, 3912, 5654, 2334, 2343, 3972, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 5726, 0, 3, 3915, 5657, 2343, 2358, 3993, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_52(buffer, 5753, 0, 3, 3924, 5669, 2358, 2373, 4014, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_53(buffer, 5780, 0, 3, 3933, 5681, 2373, 2394, 4041, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_54(buffer, 5813, 0, 3, 3942, 5693, 2394, 2415, 4068, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_50(buffer, 5840, 0, 3, 3951, 5705, 2415, 2436, 4089, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_39(buffer, 5849, 0, 3, 5651, 5654, 3972, 5726, 2460, 2478, 4125, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_40(buffer, 5873, 0, 3, 5654, 5657, 3993, 5753, 2478, 2499, 4152, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_41(buffer, 5906, 0, 3, 5657, 5669, 4014, 5780, 2499, 2538, 4209, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_42(buffer, 5981, 0, 3, 5669, 5681, 4041, 5813, 2538, 2577, 4263, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_43(buffer, 6047, 0, 3, 5681, 5693, 4068, 5840, 2577, 2616, 4314, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_24(buffer, 6071, 0, 3, 5708, 5717, 4098, 5849, 2670, 2697, 4335, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_25(buffer, 6112, 0, 3, 5717, 5726, 4125, 5873, 2697, 2724, 4369, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_26(buffer, 6153, 0, 3, 5726, 5753, 4152, 5906, 2724, 2787, 4463, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_27(buffer, 6320, 0, 3, 5753, 5780, 4209, 5981, 2787, 2850, 4563, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_28(buffer, 6462, 0, 3, 5780, 5813, 4263, 6047, 2850, 2913, 4657, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_11(buffer, 6542, 0, 3, 5849, 5873, 4369, 6153, 3018, 3096, 4835, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_12(buffer, 6710, 0, 3, 5873, 5906, 4463, 6320, 3096, 3177, 4958, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_13(buffer, 6933, 0, 3, 5906, 5981, 4563, 6462, 3177, 3258, 5084, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_2(buffer, 7116, 0, 3, 6071, 6112, 4715, 6542, 3417, 3507, 5201, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_3(buffer, 7320, 0, 3, 6112, 6153, 4835, 6710, 3507, 3597, 5339, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_4(buffer, 7524, 0, 3, 6153, 6320, 4958, 6933, 3597, 3687, 5477, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 7773, 3, 3867, 3870, 5618, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 7776, 3, 3870, 3873, 5621, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 7779, 3, 3873, 3876, 5624, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 7782, 3, 3876, 3879, 5627, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 7785, 3, 3879, 3882, 5630, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 7788, 0, 5615, 7773, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 7791, 0, 5618, 7776, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 7794, 0, 5621, 7779, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 7797, 0, 5624, 7782, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 7800, 0, 5627, 7785, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 7803, 3, 5633, 3909, 3912, 5654, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_29(buffer, 7806, 3, 5636, 3912, 3915, 5657, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 7809, 3, 5639, 3915, 3924, 5669, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 7821, 3, 5642, 3924, 3933, 5681, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 7833, 3, 5645, 3933, 3942, 5693, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 7845, 3, 5648, 3942, 3951, 5705, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_35(buffer, 7848, 0, 3, 5654, 7806, 3963, 3972, 5726, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_36(buffer, 7857, 0, 3, 5657, 7809, 3972, 3993, 5753, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_37(buffer, 7887, 0, 3, 5669, 7821, 3993, 4014, 5780, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_38(buffer, 7917, 0, 3, 5681, 7833, 4014, 4041, 5813, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_39(buffer, 7947, 0, 3, 5693, 7845, 4041, 4068, 5840, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_25(buffer, 7956, 0, 3, 7803, 7806, 5726, 7857, 4098, 4125, 5873, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_26(buffer, 7995, 0, 3, 7806, 7809, 5753, 7887, 4125, 4152, 5906, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_27(buffer, 8034, 0, 3, 7809, 7821, 5780, 7917, 4152, 4209, 5981, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_28(buffer, 8106, 0, 3, 7821, 7833, 5813, 7947, 4209, 4263, 6047, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_14(buffer, 8133, 0, 3, 7848, 7857, 5873, 7995, 4335, 4369, 6153, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_15(buffer, 8190, 0, 3, 7857, 7887, 5906, 8034, 4369, 4463, 6320, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_16(buffer, 8390, 0, 3, 7887, 7917, 5981, 8106, 4463, 4563, 6462, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_6(buffer, 8498, 0, 3, 7956, 7995, 6153, 8190, 4715, 4835, 6710, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_7(buffer, 8995, 0, 3, 7995, 8034, 6320, 8390, 4835, 4958, 6933, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_1(buffer, 9301, 0, 3, 8133, 8190, 6710, 8995, 5201, 5339, 7524, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_5(buffer, 9922, 3, 7788, 5651, 5654, 7806, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_16(buffer, 9925, 3, 7791, 5654, 5657, 7809, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 9928, 3, 7794, 5657, 5669, 7821, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 9931, 3, 7797, 5669, 5681, 7833, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_24(buffer, 9934, 3, 7800, 5681, 5693, 7845, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_5(buffer, 9937, 0, 3, 7803, 9922, 5708, 5717, 7848, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_18(buffer, 9946, 0, 3, 7806, 9925, 5717, 5726, 7857, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_19(buffer, 9955, 0, 3, 7809, 9928, 5726, 5753, 7887, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_20(buffer, 9964, 0, 3, 7821, 9931, 5753, 5780, 7917, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_21(buffer, 9973, 0, 3, 7833, 9934, 5780, 5813, 7947, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_11(buffer, 9982, 0, 3, 9922, 9925, 7857, 9955, 5849, 5873, 7995, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_12(buffer, 10009, 0, 3, 9925, 9928, 7887, 9964, 5873, 5906, 8034, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_13(buffer, 10036, 0, 3, 9928, 9931, 7917, 9973, 5906, 5981, 8106, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_5(buffer, 10063, 0, 3, 9937, 9946, 7956, 9982, 6071, 6112, 8133, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_6(buffer, 10123, 0, 3, 9946, 9955, 7995, 10009, 6112, 6153, 8190, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_7(buffer, 10183, 0, 3, 9955, 9964, 8034, 10036, 6153, 6320, 8390, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_2(buffer, 10294, 0, 3, 9982, 10009, 8190, 10183, 6542, 6710, 8995, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 10648, 0, 3, 10063, 10123, 8498, 10294, 7116, 7320, 9301, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 11944, 10648, 1296, ncols);
        }
    }

    simdtrf::transform_kk_tri(values, nvalues, buffer, 11944, nmax);
}

}  // namespace simdt2ceri
