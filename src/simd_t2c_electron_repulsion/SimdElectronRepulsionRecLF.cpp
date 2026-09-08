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
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

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

    auto buffer = CSimdMatrix(11050, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 7, 8, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 8, 9, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 9, 10, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 10, 11, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 11, 12, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 12, 13, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 13, 14, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 14, 15, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 15, 16, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 18, 21, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 21, 24, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 24, 27, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 27, 30, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 30, 33, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 33, 36, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 36, 39, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 39, 42, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 182, 0, 48, 54, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 197, 0, 54, 60, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 212, 0, 60, 66, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 227, 0, 66, 72, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 242, 0, 72, 78, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 257, 0, 78, 84, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 272, 0, 84, 90, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 287, 0, 102, 112, 197, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 308, 0, 112, 122, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 329, 0, 122, 132, 227, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 350, 0, 132, 142, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 371, 0, 142, 152, 257, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 392, 0, 152, 162, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 413, 0, 182, 197, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 441, 0, 197, 212, 329, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 469, 0, 212, 227, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 497, 0, 227, 242, 371, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 242, 257, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 553, 0, 287, 308, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 589, 0, 308, 329, 469, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 625, 0, 329, 350, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 661, 0, 350, 371, 525, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 697, 0, 413, 441, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 742, 0, 441, 469, 625, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 787, 0, 469, 497, 661, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 832, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 835, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 838, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 841, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 844, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 847, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 850, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 853, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 856, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 859, 3, 8, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 868, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 877, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 886, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 895, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 904, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 913, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 922, 3, 15, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 931, 3, 16, 45, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 940, 0, 3, 18, 859, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 958, 0, 3, 21, 868, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 976, 0, 3, 24, 877, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 994, 0, 3, 27, 886, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1012, 0, 3, 30, 895, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1030, 0, 3, 33, 904, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1048, 0, 3, 36, 913, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1066, 0, 3, 39, 922, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1084, 0, 3, 42, 931, 96, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1102, 0, 3, 54, 976, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1132, 0, 3, 60, 994, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1162, 0, 3, 66, 1012, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1192, 0, 3, 72, 1030, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1222, 0, 3, 78, 1048, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1252, 0, 3, 84, 1066, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1282, 0, 3, 90, 1084, 172, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1312, 0, 3, 102, 1102, 182, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1357, 0, 3, 112, 1132, 197, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1402, 0, 3, 122, 1162, 212, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1447, 0, 3, 132, 1192, 227, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1492, 0, 3, 142, 1222, 242, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1537, 0, 3, 152, 1252, 257, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1582, 0, 3, 162, 1282, 272, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1627, 0, 3, 197, 1402, 308, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1690, 0, 3, 212, 1447, 329, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1753, 0, 3, 227, 1492, 350, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1816, 0, 3, 242, 1537, 371, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1879, 0, 3, 257, 1582, 392, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1942, 0, 3, 287, 1627, 413, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2026, 0, 3, 308, 1690, 441, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2110, 0, 3, 329, 1753, 469, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2194, 0, 3, 350, 1816, 497, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2278, 0, 3, 371, 1879, 525, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2362, 0, 3, 441, 2110, 589, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2470, 0, 3, 469, 2194, 625, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2578, 0, 3, 497, 2278, 661, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2686, 0, 3, 553, 2362, 697, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2821, 0, 3, 589, 2470, 742, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2956, 0, 3, 625, 2578, 787, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3091, 3, 8, 9, 835, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 3097, 3, 9, 10, 838, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 3103, 3, 10, 11, 841, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3109, 3, 11, 12, 844, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3115, 3, 12, 13, 847, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3121, 3, 13, 14, 850, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3127, 3, 14, 15, 853, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3133, 3, 15, 16, 856, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3139, 0, 3, 832, 3091, 868, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3157, 0, 3, 835, 3097, 877, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3175, 0, 3, 838, 3103, 886, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3193, 0, 3, 841, 3109, 895, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3211, 0, 3, 844, 3115, 904, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3229, 0, 3, 847, 3121, 913, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3247, 0, 3, 850, 3127, 922, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3265, 0, 3, 853, 3133, 931, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3283, 0, 3, 868, 3157, 48, 54, 976,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3319, 0, 3, 877, 3175, 54, 60, 994,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3355, 0, 3, 886, 3193, 60, 66, 1012,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3391, 0, 3, 895, 3211, 66, 72, 1030,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3427, 0, 3, 904, 3229, 72, 78, 1048,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3463, 0, 3, 913, 3247, 78, 84, 1066,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3499, 0, 3, 922, 3265, 84, 90, 1084,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3535, 0, 3, 976, 3319, 102, 112, 1132,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3595, 0, 3, 994, 3355, 112, 122, 1162,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3655, 0, 3, 1012, 3391, 122, 132, 1192,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3715, 0, 3, 1030, 3427, 132, 142, 1222,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3775, 0, 3, 1048, 3463, 142, 152, 1252,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3835, 0, 3, 1066, 3499, 152, 162, 1282,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3895, 0, 3, 3283, 3319, 1132, 3595, 182,
                                                 197, 1402, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3985, 0, 3, 3319, 3355, 1162, 3655, 197,
                                                 212, 1447, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4075, 0, 3, 3355, 3391, 1192, 3715, 212,
                                                 227, 1492, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4165, 0, 3, 3391, 3427, 1222, 3775, 227,
                                                 242, 1537, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4255, 0, 3, 3427, 3463, 1252, 3835, 242,
                                                 257, 1582, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4345, 0, 3, 3535, 3595, 1402, 3985, 287,
                                                 308, 1690, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4471, 0, 3, 3595, 3655, 1447, 4075, 308,
                                                 329, 1753, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4597, 0, 3, 3655, 3715, 1492, 4165, 329,
                                                 350, 1816, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4723, 0, 3, 3715, 3775, 1537, 4255, 350,
                                                 371, 1879, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4849, 0, 3, 3895, 3985, 1690, 4471, 413,
                                                 441, 2110, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5017, 0, 3, 3985, 4075, 1753, 4597, 441,
                                                 469, 2194, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5185, 0, 3, 4075, 4165, 1816, 4723, 469,
                                                 497, 2278, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5353, 0, 3, 4345, 4471, 2110, 5017, 553,
                                                 589, 2470, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5569, 0, 3, 4471, 4597, 2194, 5185, 589,
                                                 625, 2578, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 5785, 0, 3, 4849, 5017, 2470, 5569, 697,
                                                 742, 2956, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6055, 3, 832, 835, 3097, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6065, 3, 835, 838, 3103, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6075, 3, 838, 841, 3109, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6085, 3, 841, 844, 3115, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6095, 3, 844, 847, 3121, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6105, 3, 847, 850, 3127, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6115, 3, 850, 853, 3133, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 6125, 0, 3, 3091, 6055, 3157, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6155, 0, 3, 3097, 6065, 3175, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6185, 0, 3, 3103, 6075, 3193, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6215, 0, 3, 3109, 6085, 3211, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6245, 0, 3, 3115, 6095, 3229, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6275, 0, 3, 3121, 6105, 3247, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6305, 0, 3, 3127, 6115, 3265, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 6335, 0, 3, 3139, 6125, 940, 958, 3283,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6395, 0, 3, 3157, 6155, 958, 976, 3319,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6455, 0, 3, 3175, 6185, 976, 994, 3355,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6515, 0, 3, 3193, 6215, 994, 1012, 3391,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6575, 0, 3, 3211, 6245, 1012, 1030,
                                                 3427, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6635, 0, 3, 3229, 6275, 1030, 1048,
                                                 3463, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6695, 0, 3, 3247, 6305, 1048, 1066,
                                                 3499, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6755, 0, 3, 3319, 6455, 1102, 1132,
                                                 3595, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6855, 0, 3, 3355, 6515, 1132, 1162,
                                                 3655, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6955, 0, 3, 3391, 6575, 1162, 1192,
                                                 3715, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7055, 0, 3, 3427, 6635, 1192, 1222,
                                                 3775, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7155, 0, 3, 3463, 6695, 1222, 1252,
                                                 3835, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7255, 0, 3, 6335, 6395, 3535, 6755,
                                                 1312, 1357, 3895, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7405, 0, 3, 6395, 6455, 3595, 6855,
                                                 1357, 1402, 3985, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7555, 0, 3, 6455, 6515, 3655, 6955,
                                                 1402, 1447, 4075, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7705, 0, 3, 6515, 6575, 3715, 7055,
                                                 1447, 1492, 4165, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7855, 0, 3, 6575, 6635, 3775, 7155,
                                                 1492, 1537, 4255, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8005, 0, 3, 6755, 6855, 3985, 7555,
                                                 1627, 1690, 4471, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8215, 0, 3, 6855, 6955, 4075, 7705,
                                                 1690, 1753, 4597, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8425, 0, 3, 6955, 7055, 4165, 7855,
                                                 1753, 1816, 4723, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8635, 0, 3, 7255, 7405, 4345, 8005,
                                                 1942, 2026, 4849, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8915, 0, 3, 7405, 7555, 4471, 8215,
                                                 2026, 2110, 5017, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 9195, 0, 3, 7555, 7705, 4597, 8425,
                                                 2110, 2194, 5185, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 9475, 0, 3, 8005, 8215, 5017, 9195,
                                                 2362, 2470, 5569, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 9835, 0, 3, 8635, 8915, 5353, 9475,
                                                 2686, 2821, 5785, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 10285, 9835, 450, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 10735, 10285, 45, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 10735, 7, nmax);
}

}  // namespace simdt2ceri
