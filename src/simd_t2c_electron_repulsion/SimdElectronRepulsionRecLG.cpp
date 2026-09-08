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


#include "SimdElectronRepulsionRecLG.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 20447, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 12, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 15, 16, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 16, 17, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 17, 18, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 20, 23, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 23, 26, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 26, 29, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 29, 32, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 32, 35, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 35, 38, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 38, 41, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 41, 44, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 44, 47, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 209, 0, 53, 59, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 224, 0, 59, 65, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 239, 0, 65, 71, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 254, 0, 71, 77, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 269, 0, 77, 83, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 284, 0, 83, 89, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 299, 0, 89, 95, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 314, 0, 95, 101, 189, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 329, 0, 101, 107, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 344, 0, 119, 129, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 365, 0, 129, 139, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 386, 0, 139, 149, 269, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 407, 0, 149, 159, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 428, 0, 159, 169, 299, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 449, 0, 169, 179, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 470, 0, 179, 189, 329, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 491, 0, 209, 224, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 519, 0, 224, 239, 365, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 547, 0, 239, 254, 386, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 575, 0, 254, 269, 407, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 603, 0, 269, 284, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 631, 0, 284, 299, 449, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 659, 0, 299, 314, 470, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 687, 0, 344, 365, 547, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 723, 0, 365, 386, 575, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 759, 0, 386, 407, 603, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 795, 0, 407, 428, 631, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 831, 0, 428, 449, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 867, 0, 491, 519, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 912, 0, 519, 547, 723, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 957, 0, 547, 575, 759, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1002, 0, 575, 603, 795, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1047, 0, 603, 631, 831, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1092, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1095, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1098, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1101, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1104, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1107, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1110, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1113, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1116, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1119, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1122, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1131, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1140, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1149, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1158, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1167, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1176, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1185, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1194, 3, 17, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1203, 3, 18, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1212, 0, 3, 23, 1131, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1230, 0, 3, 26, 1140, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1248, 0, 3, 29, 1149, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1266, 0, 3, 32, 1158, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1284, 0, 3, 35, 1167, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1302, 0, 3, 38, 1176, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1320, 0, 3, 41, 1185, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1338, 0, 3, 44, 1194, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1356, 0, 3, 47, 1203, 113, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1374, 0, 3, 65, 1230, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1404, 0, 3, 71, 1248, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1434, 0, 3, 77, 1266, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1464, 0, 3, 83, 1284, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1494, 0, 3, 89, 1302, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1524, 0, 3, 95, 1320, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1554, 0, 3, 101, 1338, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1584, 0, 3, 107, 1356, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1614, 0, 3, 129, 1404, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1659, 0, 3, 139, 1434, 254, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1704, 0, 3, 149, 1464, 269, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1749, 0, 3, 159, 1494, 284, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1794, 0, 3, 169, 1524, 299, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1839, 0, 3, 179, 1554, 314, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1884, 0, 3, 189, 1584, 329, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1929, 0, 3, 239, 1659, 365, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1992, 0, 3, 254, 1704, 386, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2055, 0, 3, 269, 1749, 407, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2118, 0, 3, 284, 1794, 428, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2181, 0, 3, 299, 1839, 449, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2244, 0, 3, 314, 1884, 470, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2307, 0, 3, 365, 1992, 547, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2391, 0, 3, 386, 2055, 575, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2475, 0, 3, 407, 2118, 603, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2559, 0, 3, 428, 2181, 631, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2643, 0, 3, 449, 2244, 659, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2727, 0, 3, 547, 2391, 723, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2835, 0, 3, 575, 2475, 759, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2943, 0, 3, 603, 2559, 795, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3051, 0, 3, 631, 2643, 831, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3159, 0, 3, 723, 2835, 957, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3294, 0, 3, 759, 2943, 1002, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3429, 0, 3, 795, 3051, 1047, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3564, 3, 9, 10, 1095, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3570, 3, 10, 11, 1098, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3576, 3, 11, 12, 1101, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3582, 3, 12, 13, 1104, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3588, 3, 13, 14, 1107, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3594, 3, 14, 15, 1110, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3600, 3, 15, 16, 1113, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3606, 3, 16, 17, 1116, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3612, 3, 17, 18, 1119, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3618, 0, 3, 1092, 3564, 1131, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3636, 0, 3, 1095, 3570, 1140, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3654, 0, 3, 1098, 3576, 1149, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3672, 0, 3, 1101, 3582, 1158, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3690, 0, 3, 1104, 3588, 1167, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3708, 0, 3, 1107, 3594, 1176, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3726, 0, 3, 1110, 3600, 1185, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3744, 0, 3, 1113, 3606, 1194, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3762, 0, 3, 1116, 3612, 1203, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 3780, 0, 3, 1122, 3618, 53, 59, 1212,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3816, 0, 3, 1131, 3636, 59, 65, 1230,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3852, 0, 3, 1140, 3654, 65, 71, 1248,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3888, 0, 3, 1149, 3672, 71, 77, 1266,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3924, 0, 3, 1158, 3690, 77, 83, 1284,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3960, 0, 3, 1167, 3708, 83, 89, 1302,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3996, 0, 3, 1176, 3726, 89, 95, 1320,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4032, 0, 3, 1185, 3744, 95, 101, 1338,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4068, 0, 3, 1194, 3762, 101, 107, 1356,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4104, 0, 3, 1230, 3852, 119, 129, 1404,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4164, 0, 3, 1248, 3888, 129, 139, 1434,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4224, 0, 3, 1266, 3924, 139, 149, 1464,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4284, 0, 3, 1284, 3960, 149, 159, 1494,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4344, 0, 3, 1302, 3996, 159, 169, 1524,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4404, 0, 3, 1320, 4032, 169, 179, 1554,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4464, 0, 3, 1338, 4068, 179, 189, 1584,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4524, 0, 3, 3780, 3816, 1374, 4104, 209,
                                                 224, 1614, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4614, 0, 3, 3816, 3852, 1404, 4164, 224,
                                                 239, 1659, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4704, 0, 3, 3852, 3888, 1434, 4224, 239,
                                                 254, 1704, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4794, 0, 3, 3888, 3924, 1464, 4284, 254,
                                                 269, 1749, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4884, 0, 3, 3924, 3960, 1494, 4344, 269,
                                                 284, 1794, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4974, 0, 3, 3960, 3996, 1524, 4404, 284,
                                                 299, 1839, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5064, 0, 3, 3996, 4032, 1554, 4464, 299,
                                                 314, 1884, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5154, 0, 3, 4104, 4164, 1659, 4704, 344,
                                                 365, 1992, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5280, 0, 3, 4164, 4224, 1704, 4794, 365,
                                                 386, 2055, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5406, 0, 3, 4224, 4284, 1749, 4884, 386,
                                                 407, 2118, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5532, 0, 3, 4284, 4344, 1794, 4974, 407,
                                                 428, 2181, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5658, 0, 3, 4344, 4404, 1839, 5064, 428,
                                                 449, 2244, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5784, 0, 3, 4524, 4614, 1929, 5154, 491,
                                                 519, 2307, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5952, 0, 3, 4614, 4704, 1992, 5280, 519,
                                                 547, 2391, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6120, 0, 3, 4704, 4794, 2055, 5406, 547,
                                                 575, 2475, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6288, 0, 3, 4794, 4884, 2118, 5532, 575,
                                                 603, 2559, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6456, 0, 3, 4884, 4974, 2181, 5658, 603,
                                                 631, 2643, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6624, 0, 3, 5154, 5280, 2391, 6120, 687,
                                                 723, 2835, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6840, 0, 3, 5280, 5406, 2475, 6288, 723,
                                                 759, 2943, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7056, 0, 3, 5406, 5532, 2559, 6456, 759,
                                                 795, 3051, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 7272, 0, 3, 5784, 5952, 2727, 6624, 867,
                                                 912, 3159, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 7542, 0, 3, 5952, 6120, 2835, 6840, 912,
                                                 957, 3294, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 7812, 0, 3, 6120, 6288, 2943, 7056, 957,
                                                 1002, 3429, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8082, 3, 1092, 1095, 3570, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8092, 3, 1095, 1098, 3576, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8102, 3, 1098, 1101, 3582, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8112, 3, 1101, 1104, 3588, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8122, 3, 1104, 1107, 3594, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8132, 3, 1107, 1110, 3600, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8142, 3, 1110, 1113, 3606, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8152, 3, 1113, 1116, 3612, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8162, 0, 3, 3564, 8082, 3636, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8192, 0, 3, 3570, 8092, 3654, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8222, 0, 3, 3576, 8102, 3672, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8252, 0, 3, 3582, 8112, 3690, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8282, 0, 3, 3588, 8122, 3708, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8312, 0, 3, 3594, 8132, 3726, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8342, 0, 3, 3600, 8142, 3744, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8372, 0, 3, 3606, 8152, 3762, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 8402, 0, 3, 3636, 8192, 1212, 1230,
                                                 3852, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8462, 0, 3, 3654, 8222, 1230, 1248,
                                                 3888, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8522, 0, 3, 3672, 8252, 1248, 1266,
                                                 3924, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8582, 0, 3, 3690, 8282, 1266, 1284,
                                                 3960, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8642, 0, 3, 3708, 8312, 1284, 1302,
                                                 3996, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8702, 0, 3, 3726, 8342, 1302, 1320,
                                                 4032, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8762, 0, 3, 3744, 8372, 1320, 1338,
                                                 4068, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8822, 0, 3, 3852, 8462, 1374, 1404,
                                                 4164, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8922, 0, 3, 3888, 8522, 1404, 1434,
                                                 4224, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9022, 0, 3, 3924, 8582, 1434, 1464,
                                                 4284, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9122, 0, 3, 3960, 8642, 1464, 1494,
                                                 4344, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9222, 0, 3, 3996, 8702, 1494, 1524,
                                                 4404, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9322, 0, 3, 4032, 8762, 1524, 1554,
                                                 4464, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9422, 0, 3, 8402, 8462, 4164, 8922,
                                                 1614, 1659, 4704, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9572, 0, 3, 8462, 8522, 4224, 9022,
                                                 1659, 1704, 4794, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9722, 0, 3, 8522, 8582, 4284, 9122,
                                                 1704, 1749, 4884, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9872, 0, 3, 8582, 8642, 4344, 9222,
                                                 1749, 1794, 4974, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10022, 0, 3, 8642, 8702, 4404, 9322,
                                                 1794, 1839, 5064, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10172, 0, 3, 8822, 8922, 4704, 9572,
                                                 1929, 1992, 5280, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10382, 0, 3, 8922, 9022, 4794, 9722,
                                                 1992, 2055, 5406, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10592, 0, 3, 9022, 9122, 4884, 9872,
                                                 2055, 2118, 5532, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10802, 0, 3, 9122, 9222, 4974, 10022,
                                                 2118, 2181, 5658, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11012, 0, 3, 9422, 9572, 5280, 10382,
                                                 2307, 2391, 6120, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11292, 0, 3, 9572, 9722, 5406, 10592,
                                                 2391, 2475, 6288, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11572, 0, 3, 9722, 9872, 5532, 10802,
                                                 2475, 2559, 6456, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 11852, 0, 3, 10172, 10382, 6120, 11292,
                                                 2727, 2835, 6840, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12212, 0, 3, 10382, 10592, 6288, 11572,
                                                 2835, 2943, 7056, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 12572, 0, 3, 11012, 11292, 6840, 12212,
                                                 3159, 3294, 7812, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13022, 3, 3564, 3570, 8092, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13037, 3, 3570, 3576, 8102, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13052, 3, 3576, 3582, 8112, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13067, 3, 3582, 3588, 8122, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13082, 3, 3588, 3594, 8132, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13097, 3, 3594, 3600, 8142, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13112, 3, 3600, 3606, 8152, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 13127, 0, 3, 8082, 13022, 8192, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13172, 0, 3, 8092, 13037, 8222, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13217, 0, 3, 8102, 13052, 8252, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13262, 0, 3, 8112, 13067, 8282, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13307, 0, 3, 8122, 13082, 8312, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13352, 0, 3, 8132, 13097, 8342, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13397, 0, 3, 8142, 13112, 8372, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 13442, 0, 3, 8162, 13127, 3780, 3816,
                                                 8402, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13532, 0, 3, 8192, 13172, 3816, 3852,
                                                 8462, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13622, 0, 3, 8222, 13217, 3852, 3888,
                                                 8522, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13712, 0, 3, 8252, 13262, 3888, 3924,
                                                 8582, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13802, 0, 3, 8282, 13307, 3924, 3960,
                                                 8642, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13892, 0, 3, 8312, 13352, 3960, 3996,
                                                 8702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13982, 0, 3, 8342, 13397, 3996, 4032,
                                                 8762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14072, 0, 3, 8462, 13622, 4104, 4164,
                                                 8922, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14222, 0, 3, 8522, 13712, 4164, 4224,
                                                 9022, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14372, 0, 3, 8582, 13802, 4224, 4284,
                                                 9122, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14522, 0, 3, 8642, 13892, 4284, 4344,
                                                 9222, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14672, 0, 3, 8702, 13982, 4344, 4404,
                                                 9322, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14822, 0, 3, 13442, 13532, 8822, 14072,
                                                 4524, 4614, 9422, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15047, 0, 3, 13532, 13622, 8922, 14222,
                                                 4614, 4704, 9572, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15272, 0, 3, 13622, 13712, 9022, 14372,
                                                 4704, 4794, 9722, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15497, 0, 3, 13712, 13802, 9122, 14522,
                                                 4794, 4884, 9872, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15722, 0, 3, 13802, 13892, 9222, 14672,
                                                 4884, 4974, 10022, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 15947, 0, 3, 14072, 14222, 9572, 15272,
                                                 5154, 5280, 10382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16262, 0, 3, 14222, 14372, 9722, 15497,
                                                 5280, 5406, 10592, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16577, 0, 3, 14372, 14522, 9872, 15722,
                                                 5406, 5532, 10802, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 16892, 0, 3, 14822, 15047, 10172, 15947,
                                                 5784, 5952, 11012, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 17312, 0, 3, 15047, 15272, 10382, 16262,
                                                 5952, 6120, 11292, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 17732, 0, 3, 15272, 15497, 10592, 16577,
                                                 6120, 6288, 11572, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 18152, 0, 3, 15947, 16262, 11292, 17732,
                                                 6624, 6840, 12212, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 18692, 0, 3, 16892, 17312, 11852, 18152,
                                                 7272, 7542, 12572, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 19367, 18692, 675, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 20042, 19367, 45, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 20042, 9, nmax);
}

}  // namespace simdt2ceri
