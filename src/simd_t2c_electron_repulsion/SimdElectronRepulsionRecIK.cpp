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


#include "SimdElectronRepulsionRecIK.hpp"

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
#include "SimdTransformIK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ik_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ik_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(9505, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 7, 8, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 8, 9, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 9, 10, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 10, 11, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 68, 0, 11, 12, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 71, 0, 12, 13, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 74, 0, 13, 14, 41, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 77, 0, 14, 15, 44, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 80, 0, 15, 16, 47, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 83, 0, 16, 17, 50, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 86, 0, 17, 18, 53, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 89, 0, 20, 23, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 95, 0, 23, 26, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 101, 0, 26, 29, 65, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 107, 0, 29, 32, 68, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 113, 0, 32, 35, 71, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 119, 0, 35, 38, 74, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 125, 0, 38, 41, 77, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 131, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 137, 0, 44, 47, 83, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 143, 0, 47, 50, 86, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 149, 0, 56, 59, 95, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 158, 0, 59, 62, 101, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 167, 0, 62, 65, 107, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 176, 0, 65, 68, 113, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 185, 0, 68, 71, 119, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 194, 0, 71, 74, 125, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 203, 0, 74, 77, 131, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 212, 0, 77, 80, 137, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 221, 0, 80, 83, 143, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 230, 0, 89, 95, 158, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 239, 0, 95, 101, 167, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 251, 0, 101, 107, 176, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 263, 0, 107, 113, 185, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 275, 0, 113, 119, 194, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 287, 0, 119, 125, 203, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 299, 0, 125, 131, 212, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 311, 0, 131, 137, 221, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 323, 0, 149, 158, 239, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 335, 0, 158, 167, 251, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 347, 0, 167, 176, 263, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 359, 0, 176, 185, 275, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 371, 0, 185, 194, 287, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 383, 0, 194, 203, 299, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 395, 0, 203, 212, 311, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 407, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 410, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 413, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 416, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 419, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 422, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 425, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 428, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 431, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 434, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 437, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 440, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 443, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 446, 3, 12, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 449, 3, 13, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 452, 3, 14, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 455, 3, 15, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 458, 3, 16, 47, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 461, 3, 17, 50, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 464, 3, 20, 56, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 467, 3, 23, 59, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 470, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 473, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 476, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 485, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 494, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 503, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 512, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 521, 3, 47, 83, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 524, 3, 50, 86, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 527, 3, 59, 95, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 530, 3, 62, 101, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 533, 3, 65, 107, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 542, 3, 68, 113, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 551, 3, 71, 119, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 560, 3, 74, 125, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 569, 3, 77, 131, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 578, 3, 80, 137, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 587, 3, 83, 143, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 590, 3, 89, 149, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 593, 3, 95, 158, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 596, 3, 101, 167, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 608, 3, 107, 176, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 620, 3, 113, 185, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 632, 3, 119, 194, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 644, 3, 125, 203, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 656, 3, 131, 212, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 668, 3, 137, 221, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 677, 3, 158, 239, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 692, 3, 167, 251, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 707, 3, 176, 263, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 722, 3, 185, 275, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 737, 3, 194, 287, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 752, 3, 203, 299, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 767, 3, 212, 311, ncols, p);

            compute_prim_ip_electron_repulsion_2(buffer, 782, 3, 230, 323, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 800, 3, 239, 335, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 818, 3, 251, 347, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 836, 3, 263, 359, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 854, 3, 275, 371, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 872, 3, 287, 383, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 890, 3, 299, 395, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 908, 3, 8, 9, 410, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 911, 3, 9, 10, 413, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 914, 3, 10, 11, 416, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 917, 3, 11, 12, 419, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 920, 3, 12, 13, 422, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 923, 3, 13, 14, 425, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 926, 3, 14, 15, 428, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 929, 3, 15, 16, 431, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 932, 3, 16, 17, 434, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 935, 0, 407, 908, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 938, 0, 410, 911, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 941, 0, 413, 914, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 944, 0, 416, 917, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 947, 0, 419, 920, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 950, 0, 422, 923, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 953, 0, 425, 926, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 956, 0, 428, 929, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 959, 0, 431, 932, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 962, 3, 437, 56, 59, 470, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 965, 3, 440, 59, 62, 473, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 968, 3, 443, 62, 65, 476, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 971, 3, 446, 65, 68, 485, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 980, 3, 449, 68, 71, 494, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 989, 3, 452, 71, 74, 503, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 998, 3, 455, 74, 77, 512, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1007, 3, 458, 77, 80, 521, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1010, 3, 461, 80, 83, 524, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1013, 0, 3, 470, 965, 89, 95, 530, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 1022, 0, 3, 473, 968, 95, 101, 533, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1031, 0, 3, 476, 971, 101, 107, 542, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1046, 0, 3, 485, 980, 107, 113, 551, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1061, 0, 3, 494, 989, 113, 119, 560, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1076, 0, 3, 503, 998, 119, 125, 569, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 1091, 0, 3, 512, 1007, 125, 131, 578, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1106, 0, 3, 521, 1010, 131, 137, 587, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1115, 0, 3, 962, 965, 530, 1022, 149, 158, 596, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1130, 0, 3, 965, 968, 533, 1031, 158, 167, 608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 1154, 0, 3, 968, 971, 542, 1046, 167, 176, 620, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1178, 0, 3, 971, 980, 551, 1061, 176, 185, 632, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1202, 0, 3, 980, 989, 560, 1076, 185, 194, 644, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1226, 0, 3, 989, 998, 569, 1091, 194, 203, 656, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_29(buffer, 1250, 0, 3, 998, 1007, 578, 1106, 203, 212, 668, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_15(buffer, 1271, 0, 3, 1013, 1022, 596, 1130, 230, 239, 692, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1301, 0, 3, 1022, 1031, 608, 1154, 239, 251, 707, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1334, 0, 3, 1031, 1046, 620, 1178, 251, 263, 722, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1367, 0, 3, 1046, 1061, 632, 1202, 263, 275, 737, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1400, 0, 3, 1061, 1076, 644, 1226, 275, 287, 752, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1433, 0, 3, 1076, 1091, 656, 1250, 287, 299, 767, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_5(buffer, 1466, 0, 3, 1115, 1130, 692, 1301, 323, 335, 818, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1502, 0, 3, 1130, 1154, 707, 1334, 335, 347, 836, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1538, 0, 3, 1154, 1178, 722, 1367, 347, 359, 854, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1574, 0, 3, 1178, 1202, 737, 1400, 359, 371, 872, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1610, 0, 3, 1202, 1226, 752, 1433, 371, 383, 890, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1646, 3, 407, 410, 911, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1649, 3, 410, 413, 914, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1652, 3, 413, 416, 917, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1655, 3, 416, 419, 920, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1658, 3, 419, 422, 923, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1661, 3, 422, 425, 926, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1664, 3, 425, 428, 929, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1667, 3, 428, 431, 932, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1670, 0, 911, 1649, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1673, 0, 914, 1652, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1676, 0, 917, 1655, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1679, 0, 920, 1658, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1682, 0, 923, 1661, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1685, 0, 926, 1664, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1688, 0, 929, 1667, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1691, 3, 935, 464, 467, 962, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1694, 3, 938, 467, 470, 965, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1697, 3, 941, 470, 473, 968, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_25(buffer, 1703, 3, 944, 473, 476, 971, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1709, 3, 947, 476, 485, 980, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1721, 3, 950, 485, 494, 989, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1733, 3, 953, 494, 503, 998, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_27(buffer, 1745, 3, 956, 503, 512, 1007, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_23(buffer, 1751, 3, 959, 512, 521, 1010, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_31(buffer, 1754, 0, 3, 965, 1697, 527, 530, 1022, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_32(buffer, 1763, 0, 3, 968, 1703, 530, 533, 1031, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 1775, 0, 3, 971, 1709, 533, 542, 1046, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 1802, 0, 3, 980, 1721, 542, 551, 1061, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 1829, 0, 3, 989, 1733, 551, 560, 1076, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_34(buffer, 1856, 0, 3, 998, 1745, 560, 569, 1091, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_30(buffer, 1883, 0, 3, 1007, 1751, 569, 578, 1106, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_29(buffer, 1892, 0, 3, 1691, 1694, 1013, 1754, 590, 593, 1115, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_30(buffer, 1910, 0, 3, 1694, 1697, 1022, 1763, 593, 596, 1130, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_37(buffer, 1928, 0, 3, 1697, 1703, 1031, 1775, 596, 608, 1154, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_38(buffer, 1973, 0, 3, 1703, 1709, 1046, 1802, 608, 620, 1178, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 2018, 0, 3, 1709, 1721, 1061, 1829, 620, 632, 1202, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 2063, 0, 3, 1721, 1733, 1076, 1856, 632, 644, 1226, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_33(buffer, 2108, 0, 3, 1733, 1745, 1091, 1883, 644, 656, 1250, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_20(buffer, 2141, 0, 3, 1754, 1763, 1130, 1928, 677, 692, 1301, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_21(buffer, 2201, 0, 3, 1763, 1775, 1154, 1973, 692, 707, 1334, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_22(buffer, 2264, 0, 3, 1775, 1802, 1178, 2018, 707, 722, 1367, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_22(buffer, 2327, 0, 3, 1802, 1829, 1202, 2063, 722, 737, 1400, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_19(buffer, 2390, 0, 3, 1829, 1856, 1226, 2108, 737, 752, 1433, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_2(buffer, 2450, 0, 3, 1892, 1910, 1271, 2141, 782, 800, 1466, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_5(buffer, 2522, 0, 3, 1910, 1928, 1301, 2201, 800, 818, 1502, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_8(buffer, 2594, 0, 3, 1928, 1973, 1334, 2264, 818, 836, 1538, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_8(buffer, 2666, 0, 3, 1973, 2018, 1367, 2327, 836, 854, 1574, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_9(buffer, 2738, 0, 3, 2018, 2063, 1400, 2390, 854, 872, 1610, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2810, 3, 908, 911, 1649, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2813, 3, 911, 914, 1652, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2816, 3, 914, 917, 1655, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2819, 3, 917, 920, 1658, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2822, 3, 920, 923, 1661, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2825, 3, 923, 926, 1664, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2828, 3, 926, 929, 1667, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2831, 0, 1646, 2810, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2834, 0, 1649, 2813, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2837, 0, 1652, 2816, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2840, 0, 1655, 2819, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2843, 0, 1658, 2822, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2846, 0, 1661, 2825, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2849, 0, 1664, 2828, ncols, p);

            compute_prim_dg_electron_repulsion_27(buffer, 2852, 3, 1670, 962, 965, 1697, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 2855, 3, 1673, 965, 968, 1703, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_32(buffer, 2864, 3, 1676, 968, 971, 1709, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 2873, 3, 1679, 971, 980, 1721, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 2888, 3, 1682, 980, 989, 1733, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_34(buffer, 2903, 3, 1685, 989, 998, 1745, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_26(buffer, 2912, 3, 1688, 998, 1007, 1751, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_38(buffer, 2915, 0, 3, 1697, 2855, 1013, 1022, 1763, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_39(buffer, 2930, 0, 3, 1703, 2864, 1022, 1031, 1775, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_40(buffer, 2945, 0, 3, 1709, 2873, 1031, 1046, 1802, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_41(buffer, 2987, 0, 3, 1721, 2888, 1046, 1061, 1829, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_37(buffer, 3026, 0, 3, 1733, 2903, 1061, 1076, 1856, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_33(buffer, 3065, 0, 3, 1745, 2912, 1076, 1091, 1883, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_30(buffer, 3074, 0, 3, 2852, 2855, 1763, 2930, 1115, 1130, 1928, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_34(buffer, 3095, 0, 3, 2855, 2864, 1775, 2945, 1130, 1154, 1973, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_35(buffer, 3161, 0, 3, 2864, 2873, 1802, 2987, 1154, 1178, 2018, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_36(buffer, 3233, 0, 3, 2873, 2888, 1829, 3026, 1178, 1202, 2063, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_33(buffer, 3299, 0, 3, 2888, 2903, 1856, 3065, 1202, 1226, 2108, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_16(buffer, 3344, 0, 3, 2915, 2930, 1928, 3095, 1271, 1301, 2201, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_17(buffer, 3435, 0, 3, 2930, 2945, 1973, 3161, 1301, 1334, 2264, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_18(buffer, 3529, 0, 3, 2945, 2987, 2018, 3233, 1334, 1367, 2327, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_19(buffer, 3626, 0, 3, 2987, 3026, 2063, 3299, 1367, 1400, 2390, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_5(buffer, 3714, 0, 3, 3074, 3095, 2201, 3435, 1466, 1502, 2594, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_6(buffer, 3823, 0, 3, 3095, 3161, 2264, 3529, 1502, 1538, 2666, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_7(buffer, 3932, 0, 3, 3161, 3233, 2327, 3626, 1538, 1574, 2738, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4041, 3, 1646, 1649, 2813, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4044, 3, 1649, 1652, 2816, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4047, 3, 1652, 1655, 2819, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4050, 3, 1655, 1658, 2822, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4053, 3, 1658, 1661, 2825, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4056, 3, 1661, 1664, 2828, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4059, 0, 2813, 4044, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4062, 0, 2816, 4047, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4065, 0, 2819, 4050, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4068, 0, 2822, 4053, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4071, 0, 2825, 4056, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 4074, 3, 2831, 1691, 1694, 2852, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_32(buffer, 4077, 3, 2834, 1694, 1697, 2855, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 4080, 3, 2837, 1697, 1703, 2864, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_37(buffer, 4092, 3, 2840, 1703, 1709, 2873, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_38(buffer, 4104, 3, 2843, 1709, 1721, 2888, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_39(buffer, 4122, 3, 2846, 1721, 1733, 2903, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_31(buffer, 4134, 3, 2849, 1733, 1745, 2912, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_34(buffer, 4137, 0, 3, 2855, 4080, 1754, 1763, 2930, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_35(buffer, 4146, 0, 3, 2864, 4092, 1763, 1775, 2945, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_36(buffer, 4164, 0, 3, 2873, 4104, 1775, 1802, 2987, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_37(buffer, 4221, 0, 3, 2888, 4122, 1802, 1829, 3026, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_33(buffer, 4272, 0, 3, 2903, 4134, 1829, 1856, 3065, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_24(buffer, 4281, 0, 3, 4074, 4077, 2915, 4137, 1892, 1910, 3074, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_25(buffer, 4305, 0, 3, 4077, 4080, 2930, 4146, 1910, 1928, 3095, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_26(buffer, 4329, 0, 3, 4080, 4092, 2945, 4164, 1928, 1973, 3161, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_27(buffer, 4447, 0, 3, 4092, 4104, 2987, 4221, 1973, 2018, 3233, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_28(buffer, 4547, 0, 3, 4104, 4122, 3026, 4272, 2018, 2063, 3299, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_11(buffer, 4607, 0, 3, 4137, 4146, 3095, 4329, 2141, 2201, 3435, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_12(buffer, 4732, 0, 3, 4146, 4164, 3161, 4447, 2201, 2264, 3529, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_13(buffer, 4901, 0, 3, 4164, 4221, 3233, 4547, 2264, 2327, 3626, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_2(buffer, 5035, 0, 3, 4281, 4305, 3344, 4607, 2450, 2522, 3714, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_3(buffer, 5194, 0, 3, 4305, 4329, 3435, 4732, 2522, 2594, 3823, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_4(buffer, 5353, 0, 3, 4329, 4447, 3529, 4901, 2594, 2666, 3932, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 5545, 3, 2810, 2813, 4044, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 5548, 3, 2813, 2816, 4047, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 5551, 3, 2816, 2819, 4050, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 5554, 3, 2819, 2822, 4053, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 5557, 3, 2822, 2825, 4056, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 5560, 0, 4041, 5545, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 5563, 0, 4044, 5548, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 5566, 0, 4047, 5551, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 5569, 0, 4050, 5554, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 5572, 0, 4053, 5557, ncols, p);

            compute_prim_di_electron_repulsion_29(buffer, 5575, 3, 4059, 2852, 2855, 4080, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 5578, 3, 4062, 2855, 2864, 4092, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_31(buffer, 5590, 3, 4065, 2864, 2873, 4104, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_32(buffer, 5602, 3, 4068, 2873, 2888, 4122, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_28(buffer, 5614, 3, 4071, 2888, 2903, 4134, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_23(buffer, 5617, 0, 3, 4080, 5578, 2915, 2930, 4146, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_24(buffer, 5638, 0, 3, 4092, 5590, 2930, 2945, 4164, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_25(buffer, 5659, 0, 3, 4104, 5602, 2945, 2987, 4221, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_26(buffer, 5713, 0, 3, 4122, 5614, 2987, 3026, 4272, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_14(buffer, 5722, 0, 3, 5575, 5578, 4146, 5638, 3074, 3095, 4329, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_15(buffer, 5749, 0, 3, 5578, 5590, 4164, 5659, 3095, 3161, 4447, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_16(buffer, 5895, 0, 3, 5590, 5602, 4221, 5713, 3161, 3233, 4547, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_6(buffer, 5973, 0, 3, 5617, 5638, 4329, 5749, 3344, 3435, 4732, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_7(buffer, 6345, 0, 3, 5638, 5659, 4447, 5895, 3435, 3529, 4901, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_1(buffer, 6571, 0, 3, 5722, 5749, 4732, 6345, 3714, 3823, 5353, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_5(buffer, 7051, 3, 5560, 4074, 4077, 5575, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_16(buffer, 7054, 3, 5563, 4077, 4080, 5578, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 7057, 3, 5566, 4080, 4092, 5590, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_18(buffer, 7060, 3, 5569, 4092, 4104, 5602, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_19(buffer, 7063, 3, 5572, 4104, 4122, 5614, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_11(buffer, 7066, 0, 3, 5578, 7057, 4137, 4146, 5638, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_12(buffer, 7075, 0, 3, 5590, 7060, 4146, 4164, 5659, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_13(buffer, 7084, 0, 3, 5602, 7063, 4164, 4221, 5713, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_5(buffer, 7093, 0, 3, 7051, 7054, 5617, 7066, 4281, 4305, 5722, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_6(buffer, 7120, 0, 3, 7054, 7057, 5638, 7075, 4305, 4329, 5749, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_7(buffer, 7147, 0, 3, 7057, 7060, 5659, 7084, 4329, 4447, 5895, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_2(buffer, 7225, 0, 3, 7066, 7075, 5749, 7147, 4607, 4732, 6345, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 7489, 0, 3, 7093, 7120, 5973, 7225, 5035, 5194, 6571, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8497, 7489, 1008, ncols);
        }
    }

    simdtrf::transform_ik(values, nvalues, buffer, 8497, nmax);
}

}  // namespace simdt2ceri
