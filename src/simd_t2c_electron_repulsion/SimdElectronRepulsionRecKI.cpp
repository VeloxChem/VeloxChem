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


#include "SimdElectronRepulsionRecKI.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
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
#include "SimdTransformKI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ki_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ki_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(9104, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 7, 8, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 8, 9, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 9, 10, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 68, 0, 10, 11, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 71, 0, 11, 12, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 74, 0, 12, 13, 41, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 77, 0, 13, 14, 44, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 80, 0, 14, 15, 47, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 83, 0, 15, 16, 50, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 86, 0, 16, 17, 53, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 89, 0, 17, 18, 56, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 92, 0, 20, 23, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 98, 0, 23, 26, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 104, 0, 26, 29, 65, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 110, 0, 29, 32, 68, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 116, 0, 32, 35, 71, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 122, 0, 35, 38, 74, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 128, 0, 38, 41, 77, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 134, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 140, 0, 44, 47, 83, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 146, 0, 47, 50, 86, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 152, 0, 50, 53, 89, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 158, 0, 59, 62, 104, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 167, 0, 62, 65, 110, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 176, 0, 65, 68, 116, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 185, 0, 68, 71, 122, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 194, 0, 71, 74, 128, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 203, 0, 74, 77, 134, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 212, 0, 77, 80, 140, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 221, 0, 80, 83, 146, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 230, 0, 83, 86, 152, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 239, 0, 92, 98, 158, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 251, 0, 98, 104, 167, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 263, 0, 104, 110, 176, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 275, 0, 110, 116, 185, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 287, 0, 116, 122, 194, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 299, 0, 122, 128, 203, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 311, 0, 128, 134, 212, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 323, 0, 134, 140, 221, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 335, 0, 140, 146, 230, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 347, 0, 158, 167, 263, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 362, 0, 167, 176, 275, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 377, 0, 176, 185, 287, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 392, 0, 185, 194, 299, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 407, 0, 194, 203, 311, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 422, 0, 203, 212, 323, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 437, 0, 212, 221, 335, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 452, 0, 239, 251, 347, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 467, 0, 251, 263, 362, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 482, 0, 263, 275, 377, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 497, 0, 275, 287, 392, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 512, 0, 287, 299, 407, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 527, 0, 299, 311, 422, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 542, 0, 311, 323, 437, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 557, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 560, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 563, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 566, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 569, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 572, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 575, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 578, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 581, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 584, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 587, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 590, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 593, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 596, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 599, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 602, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 605, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 608, 3, 17, 53, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 611, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 614, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 617, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 620, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 623, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 626, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 629, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 632, 3, 47, 83, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 635, 3, 50, 86, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 638, 3, 53, 89, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 641, 3, 62, 104, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 644, 3, 65, 110, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 647, 3, 68, 116, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 650, 3, 71, 122, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 659, 3, 74, 128, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 668, 3, 77, 134, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 677, 3, 80, 140, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 686, 3, 83, 146, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 689, 3, 86, 152, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 692, 3, 104, 167, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 695, 3, 110, 176, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 698, 3, 116, 185, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 710, 3, 122, 194, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 722, 3, 128, 203, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 734, 3, 134, 212, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 746, 3, 140, 221, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 755, 3, 146, 230, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 758, 3, 167, 263, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 761, 3, 176, 275, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 776, 3, 185, 287, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 791, 3, 194, 299, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 806, 3, 203, 311, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 821, 3, 212, 323, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 836, 3, 221, 335, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 845, 3, 263, 362, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 863, 3, 275, 377, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 881, 3, 287, 392, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 899, 3, 299, 407, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 917, 3, 311, 422, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 935, 3, 323, 437, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 953, 3, 362, 482, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 974, 3, 377, 497, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 995, 3, 392, 512, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1016, 3, 407, 527, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1037, 3, 422, 542, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1058, 3, 9, 10, 560, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1061, 3, 10, 11, 563, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1064, 3, 11, 12, 566, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1067, 3, 12, 13, 569, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1070, 3, 13, 14, 572, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1073, 3, 14, 15, 575, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1076, 3, 15, 16, 578, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1079, 3, 16, 17, 581, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1082, 0, 557, 1058, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1085, 0, 560, 1061, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1088, 0, 563, 1064, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1091, 0, 566, 1067, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1094, 0, 569, 1070, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1097, 0, 572, 1073, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1100, 0, 575, 1076, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1103, 0, 578, 1079, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1106, 3, 584, 59, 62, 614, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1109, 3, 587, 62, 65, 617, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1112, 3, 590, 65, 68, 620, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1115, 3, 593, 68, 71, 623, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1118, 3, 596, 71, 74, 626, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1121, 3, 599, 74, 77, 629, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1124, 3, 602, 77, 80, 632, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1127, 3, 605, 80, 83, 635, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1130, 3, 608, 83, 86, 638, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1133, 0, 3, 611, 1106, 92, 98, 641, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1142, 0, 3, 614, 1109, 98, 104, 644, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1151, 0, 3, 617, 1112, 104, 110, 647, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 1160, 0, 3, 620, 1115, 110, 116, 650, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1169, 0, 3, 623, 1118, 116, 122, 659, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1184, 0, 3, 626, 1121, 122, 128, 668, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1199, 0, 3, 629, 1124, 128, 134, 677, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1214, 0, 3, 632, 1127, 134, 140, 686, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1223, 0, 3, 635, 1130, 140, 146, 689, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1232, 0, 3, 1106, 1109, 644, 1151, 158, 167, 695, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1247, 0, 3, 1109, 1112, 647, 1160, 167, 176, 698, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1262, 0, 3, 1112, 1115, 650, 1169, 176, 185, 710, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1286, 0, 3, 1115, 1118, 659, 1184, 185, 194, 722, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1310, 0, 3, 1118, 1121, 668, 1199, 194, 203, 734, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_37(buffer, 1334, 0, 3, 1121, 1124, 677, 1214, 203, 212, 746, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1355, 0, 3, 1124, 1127, 686, 1223, 212, 221, 755, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1370, 0, 3, 1133, 1142, 692, 1232, 239, 251, 758, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1391, 0, 3, 1142, 1151, 695, 1247, 251, 263, 761, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1412, 0, 3, 1151, 1160, 698, 1262, 263, 275, 776, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1445, 0, 3, 1160, 1169, 710, 1286, 275, 287, 791, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1478, 0, 3, 1169, 1184, 722, 1310, 287, 299, 806, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1511, 0, 3, 1184, 1199, 734, 1334, 299, 311, 821, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_34(buffer, 1544, 0, 3, 1199, 1214, 746, 1355, 311, 323, 836, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 1571, 0, 3, 1232, 1247, 761, 1412, 347, 362, 863, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_20(buffer, 1613, 0, 3, 1247, 1262, 776, 1445, 362, 377, 881, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 1655, 0, 3, 1262, 1286, 791, 1478, 377, 392, 899, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 1697, 0, 3, 1286, 1310, 806, 1511, 392, 407, 917, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_18(buffer, 1739, 0, 3, 1310, 1334, 821, 1544, 407, 422, 935, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_2(buffer, 1781, 0, 3, 1370, 1391, 845, 1571, 452, 467, 953, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_5(buffer, 1826, 0, 3, 1391, 1412, 863, 1613, 467, 482, 974, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 1871, 0, 3, 1412, 1445, 881, 1655, 482, 497, 995, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 1916, 0, 3, 1445, 1478, 899, 1697, 497, 512, 1016, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 1961, 0, 3, 1478, 1511, 917, 1739, 512, 527, 1037, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2006, 3, 557, 560, 1061, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2009, 3, 560, 563, 1064, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2012, 3, 563, 566, 1067, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2015, 3, 566, 569, 1070, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2018, 3, 569, 572, 1073, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2021, 3, 572, 575, 1076, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2024, 3, 575, 578, 1079, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2027, 0, 1058, 2006, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2030, 0, 1061, 2009, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2033, 0, 1064, 2012, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2036, 0, 1067, 2015, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2039, 0, 1070, 2018, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2042, 0, 1073, 2021, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2045, 0, 1076, 2024, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2048, 3, 1082, 611, 614, 1109, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2051, 3, 1085, 614, 617, 1112, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2054, 3, 1088, 617, 620, 1115, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2060, 3, 1091, 620, 623, 1118, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2066, 3, 1094, 623, 626, 1121, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2072, 3, 1097, 626, 629, 1124, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2078, 3, 1100, 629, 632, 1127, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2084, 3, 1103, 632, 635, 1130, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2087, 0, 3, 1109, 2051, 641, 644, 1151, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2096, 0, 3, 1112, 2054, 644, 647, 1160, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_42(buffer, 2111, 0, 3, 1115, 2060, 647, 650, 1169, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2126, 0, 3, 1118, 2066, 650, 659, 1184, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2147, 0, 3, 1121, 2072, 659, 668, 1199, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_44(buffer, 2168, 0, 3, 1124, 2078, 668, 677, 1214, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_40(buffer, 2183, 0, 3, 1127, 2084, 677, 686, 1223, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_50(buffer, 2192, 0, 3, 2048, 2051, 1151, 2096, 692, 695, 1247, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_46(buffer, 2213, 0, 3, 2051, 2054, 1160, 2111, 695, 698, 1262, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_51(buffer, 2234, 0, 3, 2054, 2060, 1169, 2126, 698, 710, 1286, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2276, 0, 3, 2060, 2066, 1184, 2147, 710, 722, 1310, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_48(buffer, 2315, 0, 3, 2066, 2072, 1199, 2168, 722, 734, 1334, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_49(buffer, 2351, 0, 3, 2072, 2078, 1214, 2183, 734, 746, 1355, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_30(buffer, 2369, 0, 3, 2087, 2096, 1247, 2213, 758, 761, 1412, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_34(buffer, 2396, 0, 3, 2096, 2111, 1262, 2234, 761, 776, 1445, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_35(buffer, 2459, 0, 3, 2111, 2126, 1286, 2276, 776, 791, 1478, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_36(buffer, 2528, 0, 3, 2126, 2147, 1310, 2315, 791, 806, 1511, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_37(buffer, 2591, 0, 3, 2147, 2168, 1334, 2351, 806, 821, 1544, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_16(buffer, 2633, 0, 3, 2192, 2213, 1412, 2396, 845, 863, 1613, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_17(buffer, 2711, 0, 3, 2213, 2234, 1445, 2459, 863, 881, 1655, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_18(buffer, 2792, 0, 3, 2234, 2276, 1478, 2528, 881, 899, 1697, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_19(buffer, 2876, 0, 3, 2276, 2315, 1511, 2591, 899, 917, 1739, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_5(buffer, 2954, 0, 3, 2369, 2396, 1613, 2711, 953, 974, 1871, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_6(buffer, 3044, 0, 3, 2396, 2459, 1655, 2792, 974, 995, 1916, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_7(buffer, 3134, 0, 3, 2459, 2528, 1697, 2876, 995, 1016, 1961, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3224, 3, 1058, 1061, 2009, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3227, 3, 1061, 1064, 2012, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3230, 3, 1064, 1067, 2015, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3233, 3, 1067, 1070, 2018, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3236, 3, 1070, 1073, 2021, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3239, 3, 1073, 1076, 2024, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3242, 0, 2006, 3224, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3245, 0, 2009, 3227, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3248, 0, 2012, 3230, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3251, 0, 2015, 3233, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3254, 0, 2018, 3236, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3257, 0, 2021, 3239, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3260, 3, 2027, 1106, 1109, 2051, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 3263, 3, 2030, 1109, 1112, 2054, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3266, 3, 2033, 1112, 1115, 2060, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3275, 3, 2036, 1115, 1118, 2066, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3284, 3, 2039, 1118, 1121, 2072, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3293, 3, 2042, 1121, 1124, 2078, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3302, 3, 2045, 1124, 1127, 2084, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 3305, 0, 3, 2048, 3260, 1133, 1142, 2087, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_44(buffer, 3314, 0, 3, 2051, 3263, 1142, 1151, 2096, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 3323, 0, 3, 2054, 3266, 1151, 1160, 2111, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_50(buffer, 3344, 0, 3, 2060, 3275, 1160, 1169, 2126, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 3365, 0, 3, 2066, 3284, 1169, 1184, 2147, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_52(buffer, 3392, 0, 3, 2072, 3293, 1184, 1199, 2168, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_48(buffer, 3413, 0, 3, 2078, 3302, 1199, 1214, 2183, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_46(buffer, 3422, 0, 3, 3260, 3263, 2096, 3323, 1232, 1247, 2213, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_47(buffer, 3443, 0, 3, 3263, 3266, 2111, 3344, 1247, 1262, 2234, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_48(buffer, 3470, 0, 3, 3266, 3275, 2126, 3365, 1262, 1286, 2276, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_49(buffer, 3530, 0, 3, 3275, 3284, 2147, 3392, 1286, 1310, 2315, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_50(buffer, 3581, 0, 3, 3284, 3293, 2168, 3413, 1310, 1334, 2351, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_29(buffer, 3602, 0, 3, 3305, 3314, 2192, 3422, 1370, 1391, 2369, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_30(buffer, 3636, 0, 3, 3314, 3323, 2213, 3443, 1391, 1412, 2396, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_31(buffer, 3670, 0, 3, 3323, 3344, 2234, 3470, 1412, 1445, 2459, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_32(buffer, 3794, 0, 3, 3344, 3365, 2276, 3530, 1445, 1478, 2528, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_33(buffer, 3900, 0, 3, 3365, 3392, 2315, 3581, 1478, 1511, 2591, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_13(buffer, 3961, 0, 3, 3422, 3443, 2396, 3670, 1571, 1613, 2711, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_14(buffer, 4078, 0, 3, 3443, 3470, 2459, 3794, 1613, 1655, 2792, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_15(buffer, 4239, 0, 3, 3470, 3530, 2528, 3900, 1655, 1697, 2876, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_2(buffer, 4371, 0, 3, 3602, 3636, 2633, 3961, 1781, 1826, 2954, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_3(buffer, 4509, 0, 3, 3636, 3670, 2711, 4078, 1826, 1871, 3044, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_4(buffer, 4647, 0, 3, 3670, 3794, 2792, 4239, 1871, 1916, 3134, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4830, 3, 2006, 2009, 3227, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4833, 3, 2009, 2012, 3230, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4836, 3, 2012, 2015, 3233, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4839, 3, 2015, 2018, 3236, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4842, 3, 2018, 2021, 3239, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4845, 0, 3224, 4830, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4848, 0, 3227, 4833, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4851, 0, 3230, 4836, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4854, 0, 3233, 4839, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4857, 0, 3236, 4842, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 4860, 3, 3242, 2048, 2051, 3263, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_32(buffer, 4863, 3, 3245, 2051, 2054, 3266, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 4866, 3, 3248, 2054, 2060, 3275, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 4875, 3, 3251, 2060, 2066, 3284, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 4884, 3, 3254, 2066, 2072, 3293, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 4893, 3, 3257, 2072, 2078, 3302, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_46(buffer, 4896, 0, 3, 3263, 4863, 2087, 2096, 3323, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_47(buffer, 4905, 0, 3, 3266, 4866, 2096, 2111, 3344, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_48(buffer, 4929, 0, 3, 3275, 4875, 2111, 2126, 3365, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_49(buffer, 4953, 0, 3, 3284, 4884, 2126, 2147, 3392, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_50(buffer, 4977, 0, 3, 3293, 4893, 2147, 2168, 3413, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_35(buffer, 4986, 0, 3, 4860, 4863, 3323, 4905, 2192, 2213, 3443, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_36(buffer, 5019, 0, 3, 4863, 4866, 3344, 4929, 2213, 2234, 3470, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_37(buffer, 5052, 0, 3, 4866, 4875, 3365, 4953, 2234, 2276, 3530, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_38(buffer, 5109, 0, 3, 4875, 4884, 3392, 4977, 2276, 2315, 3581, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_21(buffer, 5133, 0, 3, 4896, 4905, 3443, 5019, 2369, 2396, 3670, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_22(buffer, 5181, 0, 3, 4905, 4929, 3470, 5052, 2396, 2459, 3794, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_23(buffer, 5334, 0, 3, 4929, 4953, 3530, 5109, 2459, 2528, 3900, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_9(buffer, 5418, 0, 3, 4986, 5019, 3670, 5181, 2633, 2711, 4078, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_10(buffer, 5792, 0, 3, 5019, 5052, 3794, 5334, 2711, 2792, 4239, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_1(buffer, 6020, 0, 3, 5133, 5181, 4078, 5792, 2954, 3044, 4647, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_9(buffer, 6497, 3, 4845, 3260, 3263, 4863, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_25(buffer, 6500, 3, 4848, 3263, 3266, 4866, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 6503, 3, 4851, 3266, 3275, 4875, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 6506, 3, 4854, 3275, 3284, 4884, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 6509, 3, 4857, 3284, 3293, 4893, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_10(buffer, 6512, 0, 3, 4860, 6497, 3305, 3314, 4896, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_31(buffer, 6521, 0, 3, 4863, 6500, 3314, 3323, 4905, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_32(buffer, 6530, 0, 3, 4866, 6503, 3323, 3344, 4929, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_33(buffer, 6539, 0, 3, 4875, 6506, 3344, 3365, 4953, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_34(buffer, 6548, 0, 3, 4884, 6509, 3365, 3392, 4977, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_22(buffer, 6557, 0, 3, 6497, 6500, 4905, 6530, 3422, 3443, 5019, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_23(buffer, 6581, 0, 3, 6500, 6503, 4929, 6539, 3443, 3470, 5052, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_24(buffer, 6605, 0, 3, 6503, 6506, 4953, 6548, 3470, 3530, 5109, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_11(buffer, 6629, 0, 3, 6512, 6521, 4986, 6557, 3602, 3636, 5133, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_12(buffer, 6680, 0, 3, 6521, 6530, 5019, 6581, 3636, 3670, 5181, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_13(buffer, 6731, 0, 3, 6530, 6539, 5052, 6605, 3670, 3794, 5334, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_5(buffer, 6818, 0, 3, 6557, 6581, 5181, 6731, 3961, 4078, 5792, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 7088, 0, 3, 6629, 6680, 5418, 6818, 4371, 4509, 6020, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8096, 7088, 1008, ncols);
        }
    }

    simdtrf::transform_ki(values, nvalues, buffer, 8096, nmax);
}

}  // namespace simdt2ceri
