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


#include "SimdElectronRepulsionRecKH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformKH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kh_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(5997, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 18, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 64, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 67, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 70, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 73, 0, 14, 15, 43, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 76, 0, 15, 16, 46, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 79, 0, 16, 17, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 19, 22, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 88, 0, 22, 25, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 94, 0, 25, 28, 61, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 28, 31, 64, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 106, 0, 31, 34, 67, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 112, 0, 34, 37, 70, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 118, 0, 37, 40, 73, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 124, 0, 40, 43, 76, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 130, 0, 43, 46, 79, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 136, 0, 52, 55, 88, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 145, 0, 55, 58, 94, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 154, 0, 58, 61, 100, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 163, 0, 61, 64, 106, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 172, 0, 64, 67, 112, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 181, 0, 67, 70, 118, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 190, 0, 70, 73, 124, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 199, 0, 73, 76, 130, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 208, 0, 82, 88, 145, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 220, 0, 88, 94, 154, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 232, 0, 94, 100, 163, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 244, 0, 100, 106, 172, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 256, 0, 106, 112, 181, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 268, 0, 112, 118, 190, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 280, 0, 118, 124, 199, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 292, 0, 136, 145, 220, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 304, 0, 145, 154, 232, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 319, 0, 154, 163, 244, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 334, 0, 163, 172, 256, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 349, 0, 172, 181, 268, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 364, 0, 181, 190, 280, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 379, 0, 208, 220, 304, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 394, 0, 220, 232, 319, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 409, 0, 232, 244, 334, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 424, 0, 244, 256, 349, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 439, 0, 256, 268, 364, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 454, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 457, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 460, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 463, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 466, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 469, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 472, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 475, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 478, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 481, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 484, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 487, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 490, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 493, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 496, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 499, 3, 16, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 502, 3, 22, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 505, 3, 25, 58, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 508, 3, 28, 61, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 511, 3, 31, 64, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 514, 3, 34, 67, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 517, 3, 37, 70, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 520, 3, 40, 73, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 523, 3, 43, 76, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 526, 3, 46, 79, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 529, 3, 52, 82, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 532, 3, 55, 88, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 535, 3, 58, 94, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 538, 3, 61, 100, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 541, 3, 64, 106, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 550, 3, 67, 112, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 559, 3, 70, 118, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 568, 3, 73, 124, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 571, 3, 76, 130, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 574, 3, 88, 145, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 577, 3, 94, 154, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 580, 3, 100, 163, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 592, 3, 106, 172, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 604, 3, 112, 181, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 616, 3, 118, 190, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 625, 3, 124, 199, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 628, 3, 136, 208, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 631, 3, 145, 220, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 634, 3, 154, 232, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 649, 3, 163, 244, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 664, 3, 172, 256, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 679, 3, 181, 268, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 694, 3, 190, 280, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 703, 3, 220, 304, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 721, 3, 232, 319, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 739, 3, 244, 334, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 757, 3, 256, 349, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 775, 3, 268, 364, ncols, p);

            compute_prim_kp_electron_repulsion_2(buffer, 793, 3, 292, 379, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 814, 3, 304, 394, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 835, 3, 319, 409, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 856, 3, 334, 424, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 877, 3, 349, 439, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 898, 3, 9, 10, 457, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 901, 3, 10, 11, 460, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 904, 3, 11, 12, 463, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 907, 3, 12, 13, 466, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 910, 3, 13, 14, 469, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 913, 3, 14, 15, 472, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 916, 3, 15, 16, 475, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 919, 0, 454, 898, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 922, 0, 457, 901, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 925, 0, 460, 904, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 928, 0, 463, 907, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 931, 0, 466, 910, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 934, 0, 469, 913, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 937, 0, 472, 916, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 940, 3, 478, 52, 55, 505, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 943, 3, 481, 55, 58, 508, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 946, 3, 484, 58, 61, 511, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 949, 3, 487, 61, 64, 514, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 952, 3, 490, 64, 67, 517, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 955, 3, 493, 67, 70, 520, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 958, 3, 496, 70, 73, 523, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 961, 3, 499, 73, 76, 526, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 964, 0, 3, 505, 943, 82, 88, 535, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 973, 0, 3, 508, 946, 88, 94, 538, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 982, 0, 3, 511, 949, 94, 100, 541, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 991, 0, 3, 514, 952, 100, 106, 550, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1006, 0, 3, 517, 955, 106, 112, 559, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1021, 0, 3, 520, 958, 112, 118, 568, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1030, 0, 3, 523, 961, 118, 124, 571, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1039, 0, 3, 940, 943, 535, 973, 136, 145, 577, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1054, 0, 3, 943, 946, 538, 982, 145, 154, 580, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1069, 0, 3, 946, 949, 541, 991, 154, 163, 592, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1093, 0, 3, 949, 952, 550, 1006, 163, 172, 604, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_37(buffer, 1117, 0, 3, 952, 955, 559, 1021, 172, 181, 616, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1138, 0, 3, 955, 958, 568, 1030, 181, 190, 625, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1153, 0, 3, 964, 973, 577, 1054, 208, 220, 634, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1174, 0, 3, 973, 982, 580, 1069, 220, 232, 649, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_33(buffer, 1207, 0, 3, 982, 991, 592, 1093, 232, 244, 664, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1243, 0, 3, 991, 1006, 604, 1117, 244, 256, 679, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_34(buffer, 1276, 0, 3, 1006, 1021, 616, 1138, 256, 268, 694, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_15(buffer, 1303, 0, 3, 1039, 1054, 634, 1174, 292, 304, 721, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_16(buffer, 1342, 0, 3, 1054, 1069, 649, 1207, 304, 319, 739, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_17(buffer, 1384, 0, 3, 1069, 1093, 664, 1243, 319, 334, 757, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_18(buffer, 1429, 0, 3, 1093, 1117, 679, 1276, 334, 349, 775, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_5(buffer, 1471, 0, 3, 1153, 1174, 721, 1342, 379, 394, 835, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_6(buffer, 1516, 0, 3, 1174, 1207, 739, 1384, 394, 409, 856, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_7(buffer, 1561, 0, 3, 1207, 1243, 757, 1429, 409, 424, 877, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1606, 3, 454, 457, 901, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1609, 3, 457, 460, 904, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1612, 3, 460, 463, 907, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1615, 3, 463, 466, 910, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1618, 3, 466, 469, 913, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1621, 3, 469, 472, 916, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1624, 0, 898, 1606, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1627, 0, 901, 1609, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1630, 0, 904, 1612, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1633, 0, 907, 1615, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1636, 0, 910, 1618, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1639, 0, 913, 1621, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1642, 3, 919, 502, 505, 943, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1645, 3, 922, 505, 508, 946, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1648, 3, 925, 508, 511, 949, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1654, 3, 928, 511, 514, 952, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1660, 3, 931, 514, 517, 955, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1666, 3, 934, 517, 520, 958, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1672, 3, 937, 520, 523, 961, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1675, 0, 3, 940, 1642, 529, 532, 964, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1684, 0, 3, 943, 1645, 532, 535, 973, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 1693, 0, 3, 946, 1648, 535, 538, 982, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_42(buffer, 1708, 0, 3, 949, 1654, 538, 541, 991, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 1723, 0, 3, 952, 1660, 541, 550, 1006, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_44(buffer, 1744, 0, 3, 955, 1666, 550, 559, 1021, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_40(buffer, 1759, 0, 3, 958, 1672, 559, 568, 1030, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_45(buffer, 1768, 0, 3, 1642, 1645, 973, 1693, 574, 577, 1054, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_46(buffer, 1786, 0, 3, 1645, 1648, 982, 1708, 577, 580, 1069, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_47(buffer, 1807, 0, 3, 1648, 1654, 991, 1723, 580, 592, 1093, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_48(buffer, 1852, 0, 3, 1654, 1660, 1006, 1744, 592, 604, 1117, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_49(buffer, 1888, 0, 3, 1660, 1666, 1021, 1759, 604, 616, 1138, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_29(buffer, 1906, 0, 3, 1675, 1684, 1039, 1768, 628, 631, 1153, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_30(buffer, 1933, 0, 3, 1684, 1693, 1054, 1786, 631, 634, 1174, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_31(buffer, 1960, 0, 3, 1693, 1708, 1069, 1807, 634, 649, 1207, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_32(buffer, 2050, 0, 3, 1708, 1723, 1093, 1852, 649, 664, 1243, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_33(buffer, 2125, 0, 3, 1723, 1744, 1117, 1888, 664, 679, 1276, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_13(buffer, 2170, 0, 3, 1768, 1786, 1174, 1960, 703, 721, 1342, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_14(buffer, 2248, 0, 3, 1786, 1807, 1207, 2050, 721, 739, 1384, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_15(buffer, 2362, 0, 3, 1807, 1852, 1243, 2125, 739, 757, 1429, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_2(buffer, 2455, 0, 3, 1906, 1933, 1303, 2170, 793, 814, 1471, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_3(buffer, 2545, 0, 3, 1933, 1960, 1342, 2248, 814, 835, 1516, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_4(buffer, 2635, 0, 3, 1960, 2050, 1384, 2362, 835, 856, 1561, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2770, 3, 898, 901, 1609, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2773, 3, 901, 904, 1612, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2776, 3, 904, 907, 1615, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2779, 3, 907, 910, 1618, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2782, 3, 910, 913, 1621, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2785, 0, 1606, 2770, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2788, 0, 1609, 2773, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2791, 0, 1612, 2776, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2794, 0, 1615, 2779, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2797, 0, 1618, 2782, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2800, 3, 1624, 940, 943, 1645, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 2803, 3, 1627, 943, 946, 1648, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 2806, 3, 1630, 946, 949, 1654, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 2812, 3, 1633, 949, 952, 1660, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 2818, 3, 1636, 952, 955, 1666, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2824, 3, 1639, 955, 958, 1672, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_44(buffer, 2827, 0, 3, 1645, 2803, 964, 973, 1693, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_45(buffer, 2836, 0, 3, 1648, 2806, 973, 982, 1708, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_46(buffer, 2854, 0, 3, 1654, 2812, 982, 991, 1723, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_47(buffer, 2872, 0, 3, 1660, 2818, 991, 1006, 1744, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_48(buffer, 2890, 0, 3, 1666, 2824, 1006, 1021, 1759, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_42(buffer, 2899, 0, 3, 2800, 2803, 1693, 2836, 1039, 1054, 1786, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_43(buffer, 2926, 0, 3, 2803, 2806, 1708, 2854, 1054, 1069, 1807, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_44(buffer, 2953, 0, 3, 2806, 2812, 1723, 2872, 1069, 1093, 1852, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_45(buffer, 2995, 0, 3, 2812, 2818, 1744, 2890, 1093, 1117, 1888, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_26(buffer, 3016, 0, 3, 2827, 2836, 1786, 2926, 1153, 1174, 1960, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_27(buffer, 3055, 0, 3, 2836, 2854, 1807, 2953, 1174, 1207, 2050, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_28(buffer, 3166, 0, 3, 2854, 2872, 1852, 2995, 1207, 1243, 2125, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_11(buffer, 3229, 0, 3, 2899, 2926, 1960, 3055, 1303, 1342, 2248, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_12(buffer, 3498, 0, 3, 2926, 2953, 2050, 3166, 1342, 1384, 2362, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_1(buffer, 3660, 0, 3, 3016, 3055, 2248, 3498, 1471, 1516, 2635, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_9(buffer, 4011, 3, 2785, 1642, 1645, 2803, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_28(buffer, 4014, 3, 2788, 1645, 1648, 2806, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 4017, 3, 2791, 1648, 1654, 2812, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 4020, 3, 2794, 1654, 1660, 2818, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 4023, 3, 2797, 1660, 1666, 2824, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_14(buffer, 4026, 0, 3, 2800, 4011, 1675, 1684, 2827, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_42(buffer, 4035, 0, 3, 2803, 4014, 1684, 1693, 2836, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_43(buffer, 4044, 0, 3, 2806, 4017, 1693, 1708, 2854, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_44(buffer, 4053, 0, 3, 2812, 4020, 1708, 1723, 2872, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_45(buffer, 4062, 0, 3, 2818, 4023, 1723, 1744, 2890, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_32(buffer, 4071, 0, 3, 4011, 4014, 2836, 4044, 1768, 1786, 2926, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_33(buffer, 4092, 0, 3, 4014, 4017, 2854, 4053, 1786, 1807, 2953, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_34(buffer, 4113, 0, 3, 4017, 4020, 2872, 4062, 1807, 1852, 2995, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_18(buffer, 4134, 0, 3, 4026, 4035, 2899, 4071, 1906, 1933, 3016, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_19(buffer, 4176, 0, 3, 4035, 4044, 2926, 4092, 1933, 1960, 3055, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_20(buffer, 4218, 0, 3, 4044, 4053, 2953, 4113, 1960, 2050, 3166, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_8(buffer, 4284, 0, 3, 4071, 4092, 3055, 4218, 2170, 2248, 3498, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 4485, 0, 3, 4134, 4176, 3229, 4284, 2455, 2545, 3660, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5241, 4485, 756, ncols);
        }
    }

    simdtrf::transform_kh(values, nvalues, buffer, 5241, nmax);
}

}  // namespace simdt2ceri
