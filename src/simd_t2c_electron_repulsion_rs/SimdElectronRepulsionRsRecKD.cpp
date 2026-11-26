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


#include "SimdElectronRepulsionRsRecKD.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_kd_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_kd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 6874, 6262, 432, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 16, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 7, 8, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 8, 9, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 9, 10, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 10, 11, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 11, 12, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 12, 13, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 13, 14, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 17, 18, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 18, 19, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 19, 20, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 20, 21, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 21, 22, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 22, 23, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 44, 47, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 53, 56, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 56, 59, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 59, 62, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 62, 65, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 65, 68, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 68, 71, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 71, 74, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 304, 0, 80, 86, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 319, 0, 86, 92, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 334, 0, 92, 98, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 349, 0, 98, 104, 214, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 364, 0, 104, 110, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 379, 0, 122, 128, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 394, 0, 128, 134, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 409, 0, 134, 140, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 424, 0, 140, 146, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 439, 0, 146, 152, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 454, 0, 164, 174, 304, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 475, 0, 174, 184, 319, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 496, 0, 184, 194, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 517, 0, 194, 204, 349, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 538, 0, 204, 214, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 559, 0, 234, 244, 379, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 580, 0, 244, 254, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 601, 0, 254, 264, 409, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 622, 0, 264, 274, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 643, 0, 274, 284, 439, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 664, 0, 304, 319, 496, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 692, 0, 319, 334, 517, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 720, 0, 334, 349, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 748, 0, 379, 394, 601, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 776, 0, 394, 409, 622, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 804, 0, 409, 424, 643, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 832, 0, 454, 475, 664, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 868, 0, 475, 496, 692, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 904, 0, 496, 517, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 940, 0, 559, 580, 748, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 976, 0, 580, 601, 776, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1012, 0, 601, 622, 804, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1048, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1051, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1054, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1057, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1060, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1063, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1066, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1069, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1072, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1075, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1078, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1081, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1084, 3, 9, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1093, 3, 10, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1102, 3, 11, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1111, 3, 12, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1120, 3, 13, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1129, 3, 14, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1138, 3, 19, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1147, 3, 20, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1156, 3, 21, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1165, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1174, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1183, 3, 24, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1192, 0, 3, 32, 1084, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1210, 0, 3, 35, 1093, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1228, 0, 3, 38, 1102, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1246, 0, 3, 41, 1111, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1264, 0, 3, 44, 1120, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1282, 0, 3, 47, 1129, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1300, 0, 3, 59, 1138, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1318, 0, 3, 62, 1147, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1336, 0, 3, 65, 1156, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1354, 0, 3, 68, 1165, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1372, 0, 3, 71, 1174, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1390, 0, 3, 74, 1183, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1408, 0, 3, 86, 1210, 184, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1438, 0, 3, 92, 1228, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1468, 0, 3, 98, 1246, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1498, 0, 3, 104, 1264, 214, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1528, 0, 3, 110, 1282, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1558, 0, 3, 128, 1318, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1588, 0, 3, 134, 1336, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1618, 0, 3, 140, 1354, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1648, 0, 3, 146, 1372, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1678, 0, 3, 152, 1390, 294, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1708, 0, 3, 184, 1438, 319, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1753, 0, 3, 194, 1468, 334, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1798, 0, 3, 204, 1498, 349, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1843, 0, 3, 214, 1528, 364, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1888, 0, 3, 254, 1588, 394, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1933, 0, 3, 264, 1618, 409, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1978, 0, 3, 274, 1648, 424, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2023, 0, 3, 284, 1678, 439, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2068, 0, 3, 319, 1753, 496, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2131, 0, 3, 334, 1798, 517, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2194, 0, 3, 349, 1843, 538, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2257, 0, 3, 394, 1933, 601, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2320, 0, 3, 409, 1978, 622, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2383, 0, 3, 424, 2023, 643, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2446, 0, 3, 496, 2131, 692, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2530, 0, 3, 517, 2194, 720, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2614, 0, 3, 601, 2320, 776, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2698, 0, 3, 622, 2383, 804, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2782, 0, 3, 692, 2530, 904, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2890, 0, 3, 776, 2698, 1012, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2998, 3, 9, 10, 1051, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3004, 3, 10, 11, 1054, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3010, 3, 11, 12, 1057, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3016, 3, 12, 13, 1060, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3022, 3, 13, 14, 1063, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3028, 3, 19, 20, 1069, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3034, 3, 20, 21, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3040, 3, 21, 22, 1075, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3046, 3, 22, 23, 1078, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3052, 3, 23, 24, 1081, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3058, 0, 3, 1048, 2998, 1093, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3076, 0, 3, 1051, 3004, 1102, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3094, 0, 3, 1054, 3010, 1111, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3112, 0, 3, 1057, 3016, 1120, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3130, 0, 3, 1060, 3022, 1129, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3148, 0, 3, 1066, 3028, 1147, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3166, 0, 3, 1069, 3034, 1156, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3184, 0, 3, 1072, 3040, 1165, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3202, 0, 3, 1075, 3046, 1174, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3220, 0, 3, 1078, 3052, 1183, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 3238, 0, 3, 1084, 3058, 80, 86, 1210,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3274, 0, 3, 1093, 3076, 86, 92, 1228,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3310, 0, 3, 1102, 3094, 92, 98, 1246,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3346, 0, 3, 1111, 3112, 98, 104, 1264,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3382, 0, 3, 1120, 3130, 104, 110, 1282,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3418, 0, 3, 1138, 3148, 122, 128, 1318,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3454, 0, 3, 1147, 3166, 128, 134, 1336,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3490, 0, 3, 1156, 3184, 134, 140, 1354,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3526, 0, 3, 1165, 3202, 140, 146, 1372,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3562, 0, 3, 1174, 3220, 146, 152, 1390,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3598, 0, 3, 1192, 3238, 164, 174, 1408,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3658, 0, 3, 1210, 3274, 174, 184, 1438,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3718, 0, 3, 1228, 3310, 184, 194, 1468,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3778, 0, 3, 1246, 3346, 194, 204, 1498,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3838, 0, 3, 1264, 3382, 204, 214, 1528,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3898, 0, 3, 1300, 3418, 234, 244, 1558,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3958, 0, 3, 1318, 3454, 244, 254, 1588,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4018, 0, 3, 1336, 3490, 254, 264, 1618,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4078, 0, 3, 1354, 3526, 264, 274, 1648,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4138, 0, 3, 1372, 3562, 274, 284, 1678,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4198, 0, 3, 3238, 3274, 1438, 3718, 304,
                                                 319, 1753, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4288, 0, 3, 3274, 3310, 1468, 3778, 319,
                                                 334, 1798, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4378, 0, 3, 3310, 3346, 1498, 3838, 334,
                                                 349, 1843, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4468, 0, 3, 3418, 3454, 1588, 4018, 379,
                                                 394, 1933, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4558, 0, 3, 3454, 3490, 1618, 4078, 394,
                                                 409, 1978, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4648, 0, 3, 3490, 3526, 1648, 4138, 409,
                                                 424, 2023, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4738, 0, 3, 3598, 3658, 1708, 4198, 454,
                                                 475, 2068, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4864, 0, 3, 3658, 3718, 1753, 4288, 475,
                                                 496, 2131, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4990, 0, 3, 3718, 3778, 1798, 4378, 496,
                                                 517, 2194, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5116, 0, 3, 3898, 3958, 1888, 4468, 559,
                                                 580, 2257, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5242, 0, 3, 3958, 4018, 1933, 4558, 580,
                                                 601, 2320, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5368, 0, 3, 4018, 4078, 1978, 4648, 601,
                                                 622, 2383, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5494, 0, 3, 4198, 4288, 2131, 4990, 664,
                                                 692, 2530, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5662, 0, 3, 4468, 4558, 2320, 5368, 748,
                                                 776, 2698, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5830, 0, 3, 4738, 4864, 2446, 5494, 832,
                                                 868, 2782, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6046, 0, 3, 5116, 5242, 2614, 5662, 940,
                                                 976, 2890, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 6262, 5830, 432, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 6694, 6478, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 6694, 5, nmax);

    simdtrf::transform_d_inner(buffer, 6694, 6262, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 75 * nvalues, nvalues, buffer, 6694, 5, nmax);
}

}  // namespace simdt2ceri
