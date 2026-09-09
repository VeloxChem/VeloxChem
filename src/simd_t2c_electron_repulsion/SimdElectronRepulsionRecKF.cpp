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


#include "SimdElectronRepulsionRecKF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kf_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 7395, 6783, 360, nvalues);

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
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 17, 20, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 20, 23, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 23, 26, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 26, 29, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 29, 32, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 32, 35, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 35, 38, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 162, 0, 44, 50, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 177, 0, 50, 56, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 192, 0, 56, 62, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 207, 0, 62, 68, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 222, 0, 68, 74, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 237, 0, 74, 80, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 252, 0, 92, 102, 177, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 273, 0, 102, 112, 192, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 294, 0, 112, 122, 207, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 315, 0, 122, 132, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 336, 0, 132, 142, 237, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 357, 0, 162, 177, 273, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 385, 0, 177, 192, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 413, 0, 192, 207, 315, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 441, 0, 207, 222, 336, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 469, 0, 252, 273, 385, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 505, 0, 273, 294, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 541, 0, 294, 315, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 577, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 580, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 583, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 586, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 589, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 592, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 595, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 598, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 607, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 616, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 625, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 634, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 643, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 652, 3, 15, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 661, 0, 3, 20, 598, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 679, 0, 3, 23, 607, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 697, 0, 3, 26, 616, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 715, 0, 3, 29, 625, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 733, 0, 3, 32, 634, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 751, 0, 3, 35, 643, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 769, 0, 3, 38, 652, 86, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 787, 0, 3, 44, 661, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 817, 0, 3, 50, 679, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 847, 0, 3, 56, 697, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 877, 0, 3, 62, 715, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 907, 0, 3, 68, 733, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 937, 0, 3, 74, 751, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 967, 0, 3, 80, 769, 152, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 997, 0, 3, 102, 847, 177, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1042, 0, 3, 112, 877, 192, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1087, 0, 3, 122, 907, 207, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1132, 0, 3, 132, 937, 222, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1177, 0, 3, 142, 967, 237, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1222, 0, 3, 162, 997, 252, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1285, 0, 3, 177, 1042, 273, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1348, 0, 3, 192, 1087, 294, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1411, 0, 3, 207, 1132, 315, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1474, 0, 3, 222, 1177, 336, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1537, 0, 3, 273, 1348, 385, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1621, 0, 3, 294, 1411, 413, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1705, 0, 3, 315, 1474, 441, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1789, 0, 3, 357, 1537, 469, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1897, 0, 3, 385, 1621, 505, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2005, 0, 3, 413, 1705, 541, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2113, 3, 9, 10, 580, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2119, 3, 10, 11, 583, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2125, 3, 11, 12, 586, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2131, 3, 12, 13, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2137, 3, 13, 14, 592, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2143, 3, 14, 15, 595, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2149, 0, 3, 577, 2113, 607, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2167, 0, 3, 580, 2119, 616, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2185, 0, 3, 583, 2125, 625, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2203, 0, 3, 586, 2131, 634, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2221, 0, 3, 589, 2137, 643, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2239, 0, 3, 592, 2143, 652, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2257, 0, 3, 598, 2149, 44, 50, 679,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2293, 0, 3, 607, 2167, 50, 56, 697,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2329, 0, 3, 616, 2185, 56, 62, 715,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2365, 0, 3, 625, 2203, 62, 68, 733,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2401, 0, 3, 634, 2221, 68, 74, 751,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2437, 0, 3, 643, 2239, 74, 80, 769,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2473, 0, 3, 679, 2293, 92, 102, 847,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2533, 0, 3, 697, 2329, 102, 112, 877,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2593, 0, 3, 715, 2365, 112, 122, 907,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2653, 0, 3, 733, 2401, 122, 132, 937,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2713, 0, 3, 751, 2437, 132, 142, 967,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2773, 0, 3, 2257, 2293, 847, 2533, 162,
                                                 177, 1042, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2863, 0, 3, 2293, 2329, 877, 2593, 177,
                                                 192, 1087, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2953, 0, 3, 2329, 2365, 907, 2653, 192,
                                                 207, 1132, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3043, 0, 3, 2365, 2401, 937, 2713, 207,
                                                 222, 1177, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3133, 0, 3, 2473, 2533, 1042, 2863, 252,
                                                 273, 1348, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3259, 0, 3, 2533, 2593, 1087, 2953, 273,
                                                 294, 1411, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3385, 0, 3, 2593, 2653, 1132, 3043, 294,
                                                 315, 1474, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3511, 0, 3, 2773, 2863, 1348, 3259, 357,
                                                 385, 1621, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3679, 0, 3, 2863, 2953, 1411, 3385, 385,
                                                 413, 1705, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 3847, 0, 3, 3133, 3259, 1621, 3679, 469,
                                                 505, 2005, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4063, 3, 577, 580, 2119, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4073, 3, 580, 583, 2125, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4083, 3, 583, 586, 2131, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4093, 3, 586, 589, 2137, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4103, 3, 589, 592, 2143, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4113, 0, 3, 2113, 4063, 2167, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4143, 0, 3, 2119, 4073, 2185, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4173, 0, 3, 2125, 4083, 2203, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4203, 0, 3, 2131, 4093, 2221, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4233, 0, 3, 2137, 4103, 2239, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4263, 0, 3, 2149, 4113, 661, 679, 2293,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4323, 0, 3, 2167, 4143, 679, 697, 2329,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4383, 0, 3, 2185, 4173, 697, 715, 2365,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4443, 0, 3, 2203, 4203, 715, 733, 2401,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4503, 0, 3, 2221, 4233, 733, 751, 2437,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4563, 0, 3, 2257, 4263, 787, 817, 2473,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4663, 0, 3, 2293, 4323, 817, 847, 2533,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4763, 0, 3, 2329, 4383, 847, 877, 2593,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4863, 0, 3, 2365, 4443, 877, 907, 2653,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4963, 0, 3, 2401, 4503, 907, 937, 2713,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5063, 0, 3, 4263, 4323, 2533, 4763, 997,
                                                 1042, 2863, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5213, 0, 3, 4323, 4383, 2593, 4863,
                                                 1042, 1087, 2953, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5363, 0, 3, 4383, 4443, 2653, 4963,
                                                 1087, 1132, 3043, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5513, 0, 3, 4563, 4663, 2773, 5063,
                                                 1222, 1285, 3133, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5723, 0, 3, 4663, 4763, 2863, 5213,
                                                 1285, 1348, 3259, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5933, 0, 3, 4763, 4863, 2953, 5363,
                                                 1348, 1411, 3385, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 6143, 0, 3, 5063, 5213, 3259, 5933,
                                                 1537, 1621, 3679, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 6423, 0, 3, 5513, 5723, 3511, 6143,
                                                 1789, 1897, 3847, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 6783, 6423, 360, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 7143, 6783, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 7143, 7, nmax);
}

}  // namespace simdt2ceri
