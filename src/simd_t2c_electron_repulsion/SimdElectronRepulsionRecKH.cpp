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
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 23322, nvalues);

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
                                            10, 11, 12}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 14, 15, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 15, 16, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 16, 17, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 19, 22, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 34, 37, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 37, 40, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 40, 43, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 43, 46, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 202, 0, 52, 58, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 217, 0, 58, 64, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 232, 0, 64, 70, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 247, 0, 70, 76, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 262, 0, 76, 82, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 277, 0, 82, 88, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 292, 0, 88, 94, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 307, 0, 94, 100, 192, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 322, 0, 112, 122, 217, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 343, 0, 122, 132, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 364, 0, 132, 142, 247, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 385, 0, 142, 152, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 406, 0, 152, 162, 277, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 427, 0, 162, 172, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 448, 0, 172, 182, 307, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 469, 0, 202, 217, 343, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 497, 0, 217, 232, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 232, 247, 385, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 553, 0, 247, 262, 406, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 262, 277, 427, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 277, 292, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 637, 0, 322, 343, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 673, 0, 343, 364, 525, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 709, 0, 364, 385, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 745, 0, 385, 406, 581, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 781, 0, 406, 427, 609, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 817, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 820, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 823, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 826, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 829, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 832, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 835, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 838, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 841, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 844, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 853, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 862, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 871, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 880, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 889, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 898, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 907, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 916, 3, 17, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 925, 0, 3, 22, 844, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 943, 0, 3, 25, 853, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 961, 0, 3, 28, 862, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 979, 0, 3, 31, 871, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 997, 0, 3, 34, 880, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1015, 0, 3, 37, 889, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1033, 0, 3, 40, 898, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1051, 0, 3, 43, 907, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1069, 0, 3, 46, 916, 106, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1087, 0, 3, 52, 925, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1117, 0, 3, 58, 943, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1147, 0, 3, 64, 961, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1177, 0, 3, 70, 979, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1207, 0, 3, 76, 997, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1237, 0, 3, 82, 1015, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1267, 0, 3, 88, 1033, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1297, 0, 3, 94, 1051, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1327, 0, 3, 100, 1069, 192, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1357, 0, 3, 122, 1147, 217, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1402, 0, 3, 132, 1177, 232, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1447, 0, 3, 142, 1207, 247, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1492, 0, 3, 152, 1237, 262, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1537, 0, 3, 162, 1267, 277, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1582, 0, 3, 172, 1297, 292, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1627, 0, 3, 182, 1327, 307, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1672, 0, 3, 202, 1357, 322, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1735, 0, 3, 217, 1402, 343, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1798, 0, 3, 232, 1447, 364, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1861, 0, 3, 247, 1492, 385, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1924, 0, 3, 262, 1537, 406, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1987, 0, 3, 277, 1582, 427, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2050, 0, 3, 292, 1627, 448, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2113, 0, 3, 343, 1798, 497, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2197, 0, 3, 364, 1861, 525, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2281, 0, 3, 385, 1924, 553, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2365, 0, 3, 406, 1987, 581, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2449, 0, 3, 427, 2050, 609, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2533, 0, 3, 469, 2113, 637, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2641, 0, 3, 497, 2197, 673, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2749, 0, 3, 525, 2281, 709, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2857, 0, 3, 553, 2365, 745, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2965, 0, 3, 581, 2449, 781, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3073, 3, 9, 10, 820, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 3079, 3, 10, 11, 823, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3085, 3, 11, 12, 826, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3091, 3, 12, 13, 829, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3097, 3, 13, 14, 832, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3103, 3, 14, 15, 835, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3109, 3, 15, 16, 838, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3115, 3, 16, 17, 841, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3121, 0, 3, 817, 3073, 853, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3139, 0, 3, 820, 3079, 862, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3157, 0, 3, 823, 3085, 871, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3175, 0, 3, 826, 3091, 880, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3193, 0, 3, 829, 3097, 889, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3211, 0, 3, 832, 3103, 898, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3229, 0, 3, 835, 3109, 907, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3247, 0, 3, 838, 3115, 916, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3265, 0, 3, 844, 3121, 52, 58, 943,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3301, 0, 3, 853, 3139, 58, 64, 961,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3337, 0, 3, 862, 3157, 64, 70, 979,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3373, 0, 3, 871, 3175, 70, 76, 997,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3409, 0, 3, 880, 3193, 76, 82, 1015,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3445, 0, 3, 889, 3211, 82, 88, 1033,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3481, 0, 3, 898, 3229, 88, 94, 1051,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3517, 0, 3, 907, 3247, 94, 100, 1069,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3553, 0, 3, 943, 3301, 112, 122, 1147,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3613, 0, 3, 961, 3337, 122, 132, 1177,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3673, 0, 3, 979, 3373, 132, 142, 1207,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3733, 0, 3, 997, 3409, 142, 152, 1237,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3793, 0, 3, 1015, 3445, 152, 162, 1267,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3853, 0, 3, 1033, 3481, 162, 172, 1297,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3913, 0, 3, 1051, 3517, 172, 182, 1327,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3973, 0, 3, 3265, 3301, 1147, 3613, 202,
                                                 217, 1402, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4063, 0, 3, 3301, 3337, 1177, 3673, 217,
                                                 232, 1447, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4153, 0, 3, 3337, 3373, 1207, 3733, 232,
                                                 247, 1492, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4243, 0, 3, 3373, 3409, 1237, 3793, 247,
                                                 262, 1537, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4333, 0, 3, 3409, 3445, 1267, 3853, 262,
                                                 277, 1582, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4423, 0, 3, 3445, 3481, 1297, 3913, 277,
                                                 292, 1627, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4513, 0, 3, 3553, 3613, 1402, 4063, 322,
                                                 343, 1798, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4639, 0, 3, 3613, 3673, 1447, 4153, 343,
                                                 364, 1861, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4765, 0, 3, 3673, 3733, 1492, 4243, 364,
                                                 385, 1924, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4891, 0, 3, 3733, 3793, 1537, 4333, 385,
                                                 406, 1987, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5017, 0, 3, 3793, 3853, 1582, 4423, 406,
                                                 427, 2050, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5143, 0, 3, 3973, 4063, 1798, 4639, 469,
                                                 497, 2197, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5311, 0, 3, 4063, 4153, 1861, 4765, 497,
                                                 525, 2281, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5479, 0, 3, 4153, 4243, 1924, 4891, 525,
                                                 553, 2365, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5647, 0, 3, 4243, 4333, 1987, 5017, 553,
                                                 581, 2449, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5815, 0, 3, 4513, 4639, 2197, 5311, 637,
                                                 673, 2749, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6031, 0, 3, 4639, 4765, 2281, 5479, 673,
                                                 709, 2857, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6247, 0, 3, 4765, 4891, 2365, 5647, 709,
                                                 745, 2965, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6463, 3, 817, 820, 3079, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6473, 3, 820, 823, 3085, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6483, 3, 823, 826, 3091, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6493, 3, 826, 829, 3097, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6503, 3, 829, 832, 3103, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6513, 3, 832, 835, 3109, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6523, 3, 835, 838, 3115, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 6533, 0, 3, 3073, 6463, 3139, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6563, 0, 3, 3079, 6473, 3157, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6593, 0, 3, 3085, 6483, 3175, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6623, 0, 3, 3091, 6493, 3193, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6653, 0, 3, 3097, 6503, 3211, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6683, 0, 3, 3103, 6513, 3229, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6713, 0, 3, 3109, 6523, 3247, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 6743, 0, 3, 3121, 6533, 925, 943, 3301,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6803, 0, 3, 3139, 6563, 943, 961, 3337,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6863, 0, 3, 3157, 6593, 961, 979, 3373,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6923, 0, 3, 3175, 6623, 979, 997, 3409,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6983, 0, 3, 3193, 6653, 997, 1015, 3445,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7043, 0, 3, 3211, 6683, 1015, 1033,
                                                 3481, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7103, 0, 3, 3229, 6713, 1033, 1051,
                                                 3517, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7163, 0, 3, 3265, 6743, 1087, 1117,
                                                 3553, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7263, 0, 3, 3301, 6803, 1117, 1147,
                                                 3613, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7363, 0, 3, 3337, 6863, 1147, 1177,
                                                 3673, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7463, 0, 3, 3373, 6923, 1177, 1207,
                                                 3733, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7563, 0, 3, 3409, 6983, 1207, 1237,
                                                 3793, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7663, 0, 3, 3445, 7043, 1237, 1267,
                                                 3853, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7763, 0, 3, 3481, 7103, 1267, 1297,
                                                 3913, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7863, 0, 3, 6743, 6803, 3613, 7363,
                                                 1357, 1402, 4063, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8013, 0, 3, 6803, 6863, 3673, 7463,
                                                 1402, 1447, 4153, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8163, 0, 3, 6863, 6923, 3733, 7563,
                                                 1447, 1492, 4243, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8313, 0, 3, 6923, 6983, 3793, 7663,
                                                 1492, 1537, 4333, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8463, 0, 3, 6983, 7043, 3853, 7763,
                                                 1537, 1582, 4423, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8613, 0, 3, 7163, 7263, 3973, 7863,
                                                 1672, 1735, 4513, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8823, 0, 3, 7263, 7363, 4063, 8013,
                                                 1735, 1798, 4639, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9033, 0, 3, 7363, 7463, 4153, 8163,
                                                 1798, 1861, 4765, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9243, 0, 3, 7463, 7563, 4243, 8313,
                                                 1861, 1924, 4891, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9453, 0, 3, 7563, 7663, 4333, 8463,
                                                 1924, 1987, 5017, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 9663, 0, 3, 7863, 8013, 4639, 9033,
                                                 2113, 2197, 5311, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 9943, 0, 3, 8013, 8163, 4765, 9243,
                                                 2197, 2281, 5479, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 10223, 0, 3, 8163, 8313, 4891, 9453,
                                                 2281, 2365, 5647, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 10503, 0, 3, 8613, 8823, 5143, 9663,
                                                 2533, 2641, 5815, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 10863, 0, 3, 8823, 9033, 5311, 9943,
                                                 2641, 2749, 6031, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 11223, 0, 3, 9033, 9243, 5479, 10223,
                                                 2749, 2857, 6247, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11583, 3, 3073, 3079, 6473, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11598, 3, 3079, 3085, 6483, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11613, 3, 3085, 3091, 6493, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11628, 3, 3091, 3097, 6503, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11643, 3, 3097, 3103, 6513, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11658, 3, 3103, 3109, 6523, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 11673, 0, 3, 6463, 11583, 6563, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11718, 0, 3, 6473, 11598, 6593, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11763, 0, 3, 6483, 11613, 6623, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11808, 0, 3, 6493, 11628, 6653, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11853, 0, 3, 6503, 11643, 6683, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11898, 0, 3, 6513, 11658, 6713, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 11943, 0, 3, 6533, 11673, 3265, 3301,
                                                 6803, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12033, 0, 3, 6563, 11718, 3301, 3337,
                                                 6863, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12123, 0, 3, 6593, 11763, 3337, 3373,
                                                 6923, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12213, 0, 3, 6623, 11808, 3373, 3409,
                                                 6983, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12303, 0, 3, 6653, 11853, 3409, 3445,
                                                 7043, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12393, 0, 3, 6683, 11898, 3445, 3481,
                                                 7103, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12483, 0, 3, 6803, 12033, 3553, 3613,
                                                 7363, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12633, 0, 3, 6863, 12123, 3613, 3673,
                                                 7463, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12783, 0, 3, 6923, 12213, 3673, 3733,
                                                 7563, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12933, 0, 3, 6983, 12303, 3733, 3793,
                                                 7663, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13083, 0, 3, 7043, 12393, 3793, 3853,
                                                 7763, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13233, 0, 3, 11943, 12033, 7363, 12633,
                                                 3973, 4063, 8013, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13458, 0, 3, 12033, 12123, 7463, 12783,
                                                 4063, 4153, 8163, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13683, 0, 3, 12123, 12213, 7563, 12933,
                                                 4153, 4243, 8313, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13908, 0, 3, 12213, 12303, 7663, 13083,
                                                 4243, 4333, 8463, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14133, 0, 3, 12483, 12633, 8013, 13458,
                                                 4513, 4639, 9033, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14448, 0, 3, 12633, 12783, 8163, 13683,
                                                 4639, 4765, 9243, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14763, 0, 3, 12783, 12933, 8313, 13908,
                                                 4765, 4891, 9453, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 15078, 0, 3, 13233, 13458, 9033, 14448,
                                                 5143, 5311, 9943, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 15498, 0, 3, 13458, 13683, 9243, 14763,
                                                 5311, 5479, 10223, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 15918, 0, 3, 14133, 14448, 9943, 15498,
                                                 5815, 6031, 11223, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16458, 3, 6463, 6473, 11598, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16479, 3, 6473, 6483, 11613, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16500, 3, 6483, 6493, 11628, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16521, 3, 6493, 6503, 11643, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16542, 3, 6503, 6513, 11658, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 16563, 0, 3, 11583, 16458, 11718, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16626, 0, 3, 11598, 16479, 11763, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16689, 0, 3, 11613, 16500, 11808, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16752, 0, 3, 11628, 16521, 11853, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16815, 0, 3, 11643, 16542, 11898, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 16878, 0, 3, 11673, 16563, 6743, 6803,
                                                 12033, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17004, 0, 3, 11718, 16626, 6803, 6863,
                                                 12123, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17130, 0, 3, 11763, 16689, 6863, 6923,
                                                 12213, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17256, 0, 3, 11808, 16752, 6923, 6983,
                                                 12303, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17382, 0, 3, 11853, 16815, 6983, 7043,
                                                 12393, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 17508, 0, 3, 11943, 16878, 7163, 7263,
                                                 12483, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 17718, 0, 3, 12033, 17004, 7263, 7363,
                                                 12633, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 17928, 0, 3, 12123, 17130, 7363, 7463,
                                                 12783, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18138, 0, 3, 12213, 17256, 7463, 7563,
                                                 12933, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18348, 0, 3, 12303, 17382, 7563, 7663,
                                                 13083, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 18558, 0, 3, 16878, 17004, 12633, 17928,
                                                 7863, 8013, 13458, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 18873, 0, 3, 17004, 17130, 12783, 18138,
                                                 8013, 8163, 13683, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19188, 0, 3, 17130, 17256, 12933, 18348,
                                                 8163, 8313, 13908, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 19503, 0, 3, 17508, 17718, 13233, 18558,
                                                 8613, 8823, 14133, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 19944, 0, 3, 17718, 17928, 13458, 18873,
                                                 8823, 9033, 14448, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 20385, 0, 3, 17928, 18138, 13683, 19188,
                                                 9033, 9243, 14763, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 20826, 0, 3, 18558, 18873, 14448, 20385,
                                                 9663, 9943, 15498, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 21414, 0, 3, 19503, 19944, 15078, 20826,
                                                 10503, 10863, 15918, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 22170, 21414, 756, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 22926, 22170, 36, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 22926, 11, nmax);
}

}  // namespace simdt2ceri
