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


#include "SimdElectronRepulsionGeom10RecLD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_ld_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_ld_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 8714, 7679, 810, nvalues);

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

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 18, 21, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 21, 24, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 24, 27, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 27, 30, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 30, 33, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 33, 36, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 36, 39, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 39, 42, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 42, 45, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 195, 0, 51, 57, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 210, 0, 57, 63, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 225, 0, 63, 69, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 240, 0, 69, 75, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 255, 0, 75, 81, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 270, 0, 81, 87, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 285, 0, 87, 93, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 300, 0, 105, 115, 195, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 321, 0, 115, 125, 210, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 342, 0, 125, 135, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 363, 0, 135, 145, 240, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 384, 0, 145, 155, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 405, 0, 155, 165, 270, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 426, 0, 165, 175, 285, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 447, 0, 195, 210, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 475, 0, 210, 225, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 503, 0, 225, 240, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 531, 0, 240, 255, 405, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 559, 0, 255, 270, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 587, 0, 300, 321, 447, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 623, 0, 321, 342, 475, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 659, 0, 342, 363, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 695, 0, 363, 384, 531, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 731, 0, 384, 405, 559, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 767, 0, 447, 475, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 812, 0, 475, 503, 695, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 857, 0, 503, 531, 731, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 902, 0, 587, 623, 767, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 957, 0, 623, 659, 812, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1012, 0, 659, 695, 857, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1067, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1070, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1073, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1076, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1079, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1082, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1085, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1088, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1091, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1100, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1109, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1118, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1127, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1136, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1145, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1154, 3, 16, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1163, 0, 3, 24, 1091, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1181, 0, 3, 27, 1100, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1199, 0, 3, 30, 1109, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1217, 0, 3, 33, 1118, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1235, 0, 3, 36, 1127, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1253, 0, 3, 39, 1136, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1271, 0, 3, 42, 1145, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1289, 0, 3, 45, 1154, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1307, 0, 3, 57, 1181, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1337, 0, 3, 63, 1199, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1367, 0, 3, 69, 1217, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1397, 0, 3, 75, 1235, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1427, 0, 3, 81, 1253, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1457, 0, 3, 87, 1271, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1487, 0, 3, 93, 1289, 185, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1517, 0, 3, 125, 1337, 210, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1562, 0, 3, 135, 1367, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1607, 0, 3, 145, 1397, 240, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1652, 0, 3, 155, 1427, 255, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1697, 0, 3, 165, 1457, 270, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1742, 0, 3, 175, 1487, 285, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1787, 0, 3, 210, 1562, 342, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1850, 0, 3, 225, 1607, 363, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1913, 0, 3, 240, 1652, 384, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1976, 0, 3, 255, 1697, 405, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2039, 0, 3, 270, 1742, 426, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2102, 0, 3, 342, 1850, 475, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2186, 0, 3, 363, 1913, 503, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2270, 0, 3, 384, 1976, 531, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2354, 0, 3, 405, 2039, 559, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2438, 0, 3, 475, 2186, 659, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2546, 0, 3, 503, 2270, 695, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2654, 0, 3, 531, 2354, 731, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2762, 0, 3, 659, 2546, 812, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2897, 0, 3, 695, 2654, 857, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 3032, 0, 3, 812, 2897, 1012, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3197, 3, 9, 10, 1070, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3203, 3, 10, 11, 1073, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3209, 3, 11, 12, 1076, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3215, 3, 12, 13, 1079, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3221, 3, 13, 14, 1082, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3227, 3, 14, 15, 1085, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3233, 3, 15, 16, 1088, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3239, 0, 3, 1067, 3197, 1100, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3257, 0, 3, 1070, 3203, 1109, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3275, 0, 3, 1073, 3209, 1118, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3293, 0, 3, 1076, 3215, 1127, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3311, 0, 3, 1079, 3221, 1136, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3329, 0, 3, 1082, 3227, 1145, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3347, 0, 3, 1085, 3233, 1154, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 3365, 0, 3, 1091, 3239, 51, 57, 1181,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3401, 0, 3, 1100, 3257, 57, 63, 1199,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3437, 0, 3, 1109, 3275, 63, 69, 1217,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3473, 0, 3, 1118, 3293, 69, 75, 1235,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3509, 0, 3, 1127, 3311, 75, 81, 1253,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3545, 0, 3, 1136, 3329, 81, 87, 1271,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3581, 0, 3, 1145, 3347, 87, 93, 1289,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3617, 0, 3, 1163, 3365, 105, 115, 1307,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3677, 0, 3, 1181, 3401, 115, 125, 1337,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3737, 0, 3, 1199, 3437, 125, 135, 1367,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3797, 0, 3, 1217, 3473, 135, 145, 1397,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3857, 0, 3, 1235, 3509, 145, 155, 1427,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3917, 0, 3, 1253, 3545, 155, 165, 1457,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3977, 0, 3, 1271, 3581, 165, 175, 1487,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4037, 0, 3, 3365, 3401, 1337, 3737, 195,
                                                 210, 1562, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4127, 0, 3, 3401, 3437, 1367, 3797, 210,
                                                 225, 1607, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4217, 0, 3, 3437, 3473, 1397, 3857, 225,
                                                 240, 1652, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4307, 0, 3, 3473, 3509, 1427, 3917, 240,
                                                 255, 1697, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4397, 0, 3, 3509, 3545, 1457, 3977, 255,
                                                 270, 1742, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4487, 0, 3, 3617, 3677, 1517, 4037, 300,
                                                 321, 1787, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4613, 0, 3, 3677, 3737, 1562, 4127, 321,
                                                 342, 1850, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4739, 0, 3, 3737, 3797, 1607, 4217, 342,
                                                 363, 1913, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4865, 0, 3, 3797, 3857, 1652, 4307, 363,
                                                 384, 1976, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4991, 0, 3, 3857, 3917, 1697, 4397, 384,
                                                 405, 2039, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5117, 0, 3, 4037, 4127, 1850, 4739, 447,
                                                 475, 2186, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5285, 0, 3, 4127, 4217, 1913, 4865, 475,
                                                 503, 2270, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5453, 0, 3, 4217, 4307, 1976, 4991, 503,
                                                 531, 2354, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5621, 0, 3, 4487, 4613, 2102, 5117, 587,
                                                 623, 2438, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5837, 0, 3, 4613, 4739, 2186, 5285, 623,
                                                 659, 2546, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6053, 0, 3, 4739, 4865, 2270, 5453, 659,
                                                 695, 2654, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 6269, 0, 3, 5117, 5285, 2546, 6053, 767,
                                                 812, 2897, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 6539, 0, 3, 5621, 5837, 2762, 6269, 902,
                                                 957, 3032, ncols, alpha, beta, p);

            compute_prim_geom_10_ld_electron_repulsion_0(buffer, 6869, 5621, 6539, ncols,
                                                         alpha);

            compute_prim_geom_10_ld_electron_repulsion_1(buffer, 7139, 5621, 6539, ncols,
                                                         alpha);

            compute_prim_geom_10_ld_electron_repulsion_2(buffer, 7409, 5621, 6539, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 7679, 6869, 810, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 8489, 7679, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 8489, 5, nmax);

    simdtrf::transform_d_inner(buffer, 8489, 7949, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 85 * nvalues, nvalues, buffer, 8489, 5, nmax);

    simdtrf::transform_d_inner(buffer, 8489, 8219, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 170 * nvalues, nvalues, buffer, 8489, 5, nmax);
}

}  // namespace simdt2ceri
