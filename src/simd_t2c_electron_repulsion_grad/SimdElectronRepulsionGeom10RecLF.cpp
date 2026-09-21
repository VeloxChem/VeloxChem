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


#include "SimdElectronRepulsionGeom10RecLF.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_lf_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_lf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 17713, 16048, 1350, nvalues);

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

            compute_prim_ls_electron_repulsion_0(buffer, 817, 0, 469, 497, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 862, 0, 497, 525, 709, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 907, 0, 525, 553, 745, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 952, 0, 553, 581, 781, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 997, 0, 637, 673, 862, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1052, 0, 673, 709, 907, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1107, 0, 709, 745, 952, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1162, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1165, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1168, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1171, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1174, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1177, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1180, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1183, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1186, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1189, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1198, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1207, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1216, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1225, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1234, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1243, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1252, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1261, 3, 17, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1270, 0, 3, 22, 1189, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1288, 0, 3, 25, 1198, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1306, 0, 3, 28, 1207, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1324, 0, 3, 31, 1216, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1342, 0, 3, 34, 1225, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1360, 0, 3, 37, 1234, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1378, 0, 3, 40, 1243, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1396, 0, 3, 43, 1252, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1414, 0, 3, 46, 1261, 106, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1432, 0, 3, 52, 1270, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1462, 0, 3, 58, 1288, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1492, 0, 3, 64, 1306, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1522, 0, 3, 70, 1324, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1552, 0, 3, 76, 1342, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1582, 0, 3, 82, 1360, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1612, 0, 3, 88, 1378, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1642, 0, 3, 94, 1396, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1672, 0, 3, 100, 1414, 192, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1702, 0, 3, 122, 1492, 217, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1747, 0, 3, 132, 1522, 232, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1792, 0, 3, 142, 1552, 247, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1837, 0, 3, 152, 1582, 262, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1882, 0, 3, 162, 1612, 277, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1927, 0, 3, 172, 1642, 292, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1972, 0, 3, 182, 1672, 307, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2017, 0, 3, 202, 1702, 322, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2080, 0, 3, 217, 1747, 343, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2143, 0, 3, 232, 1792, 364, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2206, 0, 3, 247, 1837, 385, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2269, 0, 3, 262, 1882, 406, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2332, 0, 3, 277, 1927, 427, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2395, 0, 3, 292, 1972, 448, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2458, 0, 3, 343, 2143, 497, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2542, 0, 3, 364, 2206, 525, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2626, 0, 3, 385, 2269, 553, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2710, 0, 3, 406, 2332, 581, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2794, 0, 3, 427, 2395, 609, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2878, 0, 3, 469, 2458, 637, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2986, 0, 3, 497, 2542, 673, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3094, 0, 3, 525, 2626, 709, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3202, 0, 3, 553, 2710, 745, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3310, 0, 3, 581, 2794, 781, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3418, 0, 3, 673, 3094, 862, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3553, 0, 3, 709, 3202, 907, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3688, 0, 3, 745, 3310, 952, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 3823, 0, 3, 817, 3418, 997, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 3988, 0, 3, 862, 3553, 1052, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 4153, 0, 3, 907, 3688, 1107, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4318, 3, 9, 10, 1165, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4324, 3, 10, 11, 1168, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4330, 3, 11, 12, 1171, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4336, 3, 12, 13, 1174, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4342, 3, 13, 14, 1177, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4348, 3, 14, 15, 1180, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4354, 3, 15, 16, 1183, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4360, 3, 16, 17, 1186, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4366, 0, 3, 1162, 4318, 1198, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4384, 0, 3, 1165, 4324, 1207, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4402, 0, 3, 1168, 4330, 1216, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4420, 0, 3, 1171, 4336, 1225, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4438, 0, 3, 1174, 4342, 1234, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4456, 0, 3, 1177, 4348, 1243, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4474, 0, 3, 1180, 4354, 1252, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4492, 0, 3, 1183, 4360, 1261, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4510, 0, 3, 1189, 4366, 52, 58, 1288,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4546, 0, 3, 1198, 4384, 58, 64, 1306,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4582, 0, 3, 1207, 4402, 64, 70, 1324,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4618, 0, 3, 1216, 4420, 70, 76, 1342,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4654, 0, 3, 1225, 4438, 76, 82, 1360,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4690, 0, 3, 1234, 4456, 82, 88, 1378,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4726, 0, 3, 1243, 4474, 88, 94, 1396,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4762, 0, 3, 1252, 4492, 94, 100, 1414,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4798, 0, 3, 1288, 4546, 112, 122, 1492,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4858, 0, 3, 1306, 4582, 122, 132, 1522,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4918, 0, 3, 1324, 4618, 132, 142, 1552,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4978, 0, 3, 1342, 4654, 142, 152, 1582,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5038, 0, 3, 1360, 4690, 152, 162, 1612,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5098, 0, 3, 1378, 4726, 162, 172, 1642,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5158, 0, 3, 1396, 4762, 172, 182, 1672,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5218, 0, 3, 4510, 4546, 1492, 4858, 202,
                                                 217, 1747, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5308, 0, 3, 4546, 4582, 1522, 4918, 217,
                                                 232, 1792, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5398, 0, 3, 4582, 4618, 1552, 4978, 232,
                                                 247, 1837, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5488, 0, 3, 4618, 4654, 1582, 5038, 247,
                                                 262, 1882, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5578, 0, 3, 4654, 4690, 1612, 5098, 262,
                                                 277, 1927, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5668, 0, 3, 4690, 4726, 1642, 5158, 277,
                                                 292, 1972, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5758, 0, 3, 4798, 4858, 1747, 5308, 322,
                                                 343, 2143, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5884, 0, 3, 4858, 4918, 1792, 5398, 343,
                                                 364, 2206, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6010, 0, 3, 4918, 4978, 1837, 5488, 364,
                                                 385, 2269, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6136, 0, 3, 4978, 5038, 1882, 5578, 385,
                                                 406, 2332, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6262, 0, 3, 5038, 5098, 1927, 5668, 406,
                                                 427, 2395, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6388, 0, 3, 5218, 5308, 2143, 5884, 469,
                                                 497, 2542, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6556, 0, 3, 5308, 5398, 2206, 6010, 497,
                                                 525, 2626, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6724, 0, 3, 5398, 5488, 2269, 6136, 525,
                                                 553, 2710, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6892, 0, 3, 5488, 5578, 2332, 6262, 553,
                                                 581, 2794, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7060, 0, 3, 5758, 5884, 2542, 6556, 637,
                                                 673, 3094, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7276, 0, 3, 5884, 6010, 2626, 6724, 673,
                                                 709, 3202, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7492, 0, 3, 6010, 6136, 2710, 6892, 709,
                                                 745, 3310, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 7708, 0, 3, 6388, 6556, 3094, 7276, 817,
                                                 862, 3553, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 7978, 0, 3, 6556, 6724, 3202, 7492, 862,
                                                 907, 3688, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 8248, 0, 3, 7060, 7276, 3553, 7978, 997,
                                                 1052, 4153, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8578, 3, 1162, 1165, 4324, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8588, 3, 1165, 1168, 4330, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8598, 3, 1168, 1171, 4336, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8608, 3, 1171, 1174, 4342, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8618, 3, 1174, 1177, 4348, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8628, 3, 1177, 1180, 4354, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8638, 3, 1180, 1183, 4360, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8648, 0, 3, 4318, 8578, 4384, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8678, 0, 3, 4324, 8588, 4402, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8708, 0, 3, 4330, 8598, 4420, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8738, 0, 3, 4336, 8608, 4438, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8768, 0, 3, 4342, 8618, 4456, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8798, 0, 3, 4348, 8628, 4474, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8828, 0, 3, 4354, 8638, 4492, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 8858, 0, 3, 4366, 8648, 1270, 1288,
                                                 4546, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8918, 0, 3, 4384, 8678, 1288, 1306,
                                                 4582, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8978, 0, 3, 4402, 8708, 1306, 1324,
                                                 4618, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9038, 0, 3, 4420, 8738, 1324, 1342,
                                                 4654, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9098, 0, 3, 4438, 8768, 1342, 1360,
                                                 4690, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9158, 0, 3, 4456, 8798, 1360, 1378,
                                                 4726, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9218, 0, 3, 4474, 8828, 1378, 1396,
                                                 4762, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9278, 0, 3, 4510, 8858, 1432, 1462,
                                                 4798, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9378, 0, 3, 4546, 8918, 1462, 1492,
                                                 4858, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9478, 0, 3, 4582, 8978, 1492, 1522,
                                                 4918, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9578, 0, 3, 4618, 9038, 1522, 1552,
                                                 4978, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9678, 0, 3, 4654, 9098, 1552, 1582,
                                                 5038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9778, 0, 3, 4690, 9158, 1582, 1612,
                                                 5098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9878, 0, 3, 4726, 9218, 1612, 1642,
                                                 5158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9978, 0, 3, 8858, 8918, 4858, 9478,
                                                 1702, 1747, 5308, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10128, 0, 3, 8918, 8978, 4918, 9578,
                                                 1747, 1792, 5398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10278, 0, 3, 8978, 9038, 4978, 9678,
                                                 1792, 1837, 5488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10428, 0, 3, 9038, 9098, 5038, 9778,
                                                 1837, 1882, 5578, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10578, 0, 3, 9098, 9158, 5098, 9878,
                                                 1882, 1927, 5668, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10728, 0, 3, 9278, 9378, 5218, 9978,
                                                 2017, 2080, 5758, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10938, 0, 3, 9378, 9478, 5308, 10128,
                                                 2080, 2143, 5884, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11148, 0, 3, 9478, 9578, 5398, 10278,
                                                 2143, 2206, 6010, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11358, 0, 3, 9578, 9678, 5488, 10428,
                                                 2206, 2269, 6136, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11568, 0, 3, 9678, 9778, 5578, 10578,
                                                 2269, 2332, 6262, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11778, 0, 3, 9978, 10128, 5884, 11148,
                                                 2458, 2542, 6556, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12058, 0, 3, 10128, 10278, 6010, 11358,
                                                 2542, 2626, 6724, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12338, 0, 3, 10278, 10428, 6136, 11568,
                                                 2626, 2710, 6892, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12618, 0, 3, 10728, 10938, 6388, 11778,
                                                 2878, 2986, 7060, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12978, 0, 3, 10938, 11148, 6556, 12058,
                                                 2986, 3094, 7276, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 13338, 0, 3, 11148, 11358, 6724, 12338,
                                                 3094, 3202, 7492, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 13698, 0, 3, 11778, 12058, 7276, 13338,
                                                 3418, 3553, 7978, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 14148, 0, 3, 12618, 12978, 7708, 13698,
                                                 3823, 3988, 8248, ncols, alpha, beta, p);

            simdgeo::geom_l_x(buffer, 14698, 12618, 14148, 1, 10, ncols, alpha);

            simdgeo::geom_l_y(buffer, 15148, 12618, 14148, 1, 10, ncols, alpha);

            simdgeo::geom_l_z(buffer, 15598, 12618, 14148, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 16048, 14698, 1350, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 17398, 16048, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 17398, 7, nmax);

    simdtrf::transform_f_inner(buffer, 17398, 16498, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 119 * nvalues, nvalues, buffer, 17398, 7, nmax);

    simdtrf::transform_f_inner(buffer, 17398, 16948, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 238 * nvalues, nvalues, buffer, 17398, 7, nmax);
}

}  // namespace simdt2ceri
