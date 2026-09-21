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


#include "SimdElectronRepulsionRsRecKF.hpp"

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
compute_rs_kf_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_kf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 14532, 13560, 720, nvalues);

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
                                                9, 10}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 17, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 27, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 7, 8, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 8, 9, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 9, 10, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 10, 11, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 11, 12, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 12, 13, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 13, 14, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 14, 15, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 20, 21, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 21, 22, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 22, 23, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 23, 24, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 24, 25, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 25, 26, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 28, 31, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 31, 34, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 34, 37, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 37, 40, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 40, 43, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 43, 46, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 46, 49, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 55, 58, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 58, 61, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 61, 64, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 64, 67, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 67, 70, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 70, 73, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 73, 76, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 318, 0, 82, 88, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 333, 0, 88, 94, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 348, 0, 94, 100, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 363, 0, 100, 106, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 378, 0, 106, 112, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 393, 0, 112, 118, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 408, 0, 130, 136, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 423, 0, 136, 142, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 438, 0, 142, 148, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 453, 0, 148, 154, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 468, 0, 154, 160, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 483, 0, 160, 166, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 498, 0, 178, 188, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 519, 0, 188, 198, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 540, 0, 198, 208, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 561, 0, 208, 218, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 582, 0, 218, 228, 393, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 603, 0, 248, 258, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 624, 0, 258, 268, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 645, 0, 268, 278, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 666, 0, 278, 288, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 687, 0, 288, 298, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 708, 0, 318, 333, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 736, 0, 333, 348, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 764, 0, 348, 363, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 792, 0, 363, 378, 582, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 820, 0, 408, 423, 624, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 848, 0, 423, 438, 645, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 876, 0, 438, 453, 666, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 904, 0, 453, 468, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 932, 0, 498, 519, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 968, 0, 519, 540, 764, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1004, 0, 540, 561, 792, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1040, 0, 603, 624, 848, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1076, 0, 624, 645, 876, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1112, 0, 645, 666, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1148, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1151, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1154, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1157, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1160, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1163, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1166, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1169, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1172, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1175, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1178, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1181, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1184, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1187, 3, 27, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1190, 3, 9, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1199, 3, 10, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1208, 3, 11, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1217, 3, 12, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1226, 3, 13, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1235, 3, 14, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1244, 3, 15, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1253, 3, 20, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1262, 3, 21, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1271, 3, 22, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1280, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1289, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1298, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1307, 3, 26, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1316, 0, 3, 31, 1190, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1334, 0, 3, 34, 1199, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1352, 0, 3, 37, 1208, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1370, 0, 3, 40, 1217, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1388, 0, 3, 43, 1226, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1406, 0, 3, 46, 1235, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1424, 0, 3, 49, 1244, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1442, 0, 3, 58, 1253, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1460, 0, 3, 61, 1262, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1478, 0, 3, 64, 1271, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1496, 0, 3, 67, 1280, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1514, 0, 3, 70, 1289, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1532, 0, 3, 73, 1298, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1550, 0, 3, 76, 1307, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1568, 0, 3, 82, 1316, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1598, 0, 3, 88, 1334, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1628, 0, 3, 94, 1352, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1658, 0, 3, 100, 1370, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1688, 0, 3, 106, 1388, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1718, 0, 3, 112, 1406, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1748, 0, 3, 118, 1424, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1778, 0, 3, 130, 1442, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1808, 0, 3, 136, 1460, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1838, 0, 3, 142, 1478, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1868, 0, 3, 148, 1496, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1898, 0, 3, 154, 1514, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1928, 0, 3, 160, 1532, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1958, 0, 3, 166, 1550, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1988, 0, 3, 188, 1628, 333, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2033, 0, 3, 198, 1658, 348, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2078, 0, 3, 208, 1688, 363, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2123, 0, 3, 218, 1718, 378, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2168, 0, 3, 228, 1748, 393, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2213, 0, 3, 258, 1838, 423, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2258, 0, 3, 268, 1868, 438, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2303, 0, 3, 278, 1898, 453, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2348, 0, 3, 288, 1928, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2393, 0, 3, 298, 1958, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2438, 0, 3, 318, 1988, 498, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2501, 0, 3, 333, 2033, 519, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2564, 0, 3, 348, 2078, 540, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2627, 0, 3, 363, 2123, 561, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2690, 0, 3, 378, 2168, 582, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2753, 0, 3, 408, 2213, 603, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2816, 0, 3, 423, 2258, 624, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2879, 0, 3, 438, 2303, 645, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2942, 0, 3, 453, 2348, 666, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3005, 0, 3, 468, 2393, 687, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3068, 0, 3, 519, 2564, 736, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3152, 0, 3, 540, 2627, 764, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3236, 0, 3, 561, 2690, 792, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3320, 0, 3, 624, 2879, 848, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3404, 0, 3, 645, 2942, 876, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3488, 0, 3, 666, 3005, 904, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3572, 0, 3, 708, 3068, 932, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3680, 0, 3, 736, 3152, 968, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3788, 0, 3, 764, 3236, 1004, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3896, 0, 3, 820, 3320, 1040, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4004, 0, 3, 848, 3404, 1076, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4112, 0, 3, 876, 3488, 1112, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4220, 3, 9, 10, 1151, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4226, 3, 10, 11, 1154, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4232, 3, 11, 12, 1157, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4238, 3, 12, 13, 1160, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4244, 3, 13, 14, 1163, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4250, 3, 14, 15, 1166, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4256, 3, 20, 21, 1172, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4262, 3, 21, 22, 1175, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4268, 3, 22, 23, 1178, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4274, 3, 23, 24, 1181, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4280, 3, 24, 25, 1184, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4286, 3, 25, 26, 1187, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4292, 0, 3, 1148, 4220, 1199, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4310, 0, 3, 1151, 4226, 1208, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4328, 0, 3, 1154, 4232, 1217, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4346, 0, 3, 1157, 4238, 1226, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4364, 0, 3, 1160, 4244, 1235, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4382, 0, 3, 1163, 4250, 1244, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4400, 0, 3, 1169, 4256, 1262, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4418, 0, 3, 1172, 4262, 1271, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4436, 0, 3, 1175, 4268, 1280, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4454, 0, 3, 1178, 4274, 1289, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4472, 0, 3, 1181, 4280, 1298, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4490, 0, 3, 1184, 4286, 1307, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4508, 0, 3, 1190, 4292, 82, 88, 1334,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4544, 0, 3, 1199, 4310, 88, 94, 1352,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4580, 0, 3, 1208, 4328, 94, 100, 1370,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4616, 0, 3, 1217, 4346, 100, 106, 1388,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4652, 0, 3, 1226, 4364, 106, 112, 1406,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4688, 0, 3, 1235, 4382, 112, 118, 1424,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4724, 0, 3, 1253, 4400, 130, 136, 1460,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4760, 0, 3, 1262, 4418, 136, 142, 1478,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4796, 0, 3, 1271, 4436, 142, 148, 1496,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4832, 0, 3, 1280, 4454, 148, 154, 1514,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4868, 0, 3, 1289, 4472, 154, 160, 1532,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4904, 0, 3, 1298, 4490, 160, 166, 1550,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4940, 0, 3, 1334, 4544, 178, 188, 1628,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5000, 0, 3, 1352, 4580, 188, 198, 1658,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5060, 0, 3, 1370, 4616, 198, 208, 1688,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5120, 0, 3, 1388, 4652, 208, 218, 1718,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5180, 0, 3, 1406, 4688, 218, 228, 1748,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5240, 0, 3, 1460, 4760, 248, 258, 1838,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5300, 0, 3, 1478, 4796, 258, 268, 1868,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5360, 0, 3, 1496, 4832, 268, 278, 1898,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5420, 0, 3, 1514, 4868, 278, 288, 1928,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5480, 0, 3, 1532, 4904, 288, 298, 1958,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5540, 0, 3, 4508, 4544, 1628, 5000, 318,
                                                 333, 2033, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5630, 0, 3, 4544, 4580, 1658, 5060, 333,
                                                 348, 2078, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5720, 0, 3, 4580, 4616, 1688, 5120, 348,
                                                 363, 2123, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5810, 0, 3, 4616, 4652, 1718, 5180, 363,
                                                 378, 2168, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5900, 0, 3, 4724, 4760, 1838, 5300, 408,
                                                 423, 2258, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5990, 0, 3, 4760, 4796, 1868, 5360, 423,
                                                 438, 2303, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6080, 0, 3, 4796, 4832, 1898, 5420, 438,
                                                 453, 2348, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6170, 0, 3, 4832, 4868, 1928, 5480, 453,
                                                 468, 2393, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6260, 0, 3, 4940, 5000, 2033, 5630, 498,
                                                 519, 2564, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6386, 0, 3, 5000, 5060, 2078, 5720, 519,
                                                 540, 2627, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6512, 0, 3, 5060, 5120, 2123, 5810, 540,
                                                 561, 2690, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6638, 0, 3, 5240, 5300, 2258, 5990, 603,
                                                 624, 2879, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6764, 0, 3, 5300, 5360, 2303, 6080, 624,
                                                 645, 2942, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6890, 0, 3, 5360, 5420, 2348, 6170, 645,
                                                 666, 3005, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7016, 0, 3, 5540, 5630, 2564, 6386, 708,
                                                 736, 3152, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7184, 0, 3, 5630, 5720, 2627, 6512, 736,
                                                 764, 3236, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7352, 0, 3, 5900, 5990, 2879, 6764, 820,
                                                 848, 3404, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7520, 0, 3, 5990, 6080, 2942, 6890, 848,
                                                 876, 3488, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7688, 0, 3, 6260, 6386, 3152, 7184, 932,
                                                 968, 3788, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7904, 0, 3, 6638, 6764, 3404, 7520,
                                                 1040, 1076, 4112, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8120, 3, 1148, 1151, 4226, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8130, 3, 1151, 1154, 4232, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8140, 3, 1154, 1157, 4238, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8150, 3, 1157, 1160, 4244, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8160, 3, 1160, 1163, 4250, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8170, 3, 1169, 1172, 4262, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8180, 3, 1172, 1175, 4268, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8190, 3, 1175, 1178, 4274, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8200, 3, 1178, 1181, 4280, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8210, 3, 1181, 1184, 4286, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8220, 0, 3, 4220, 8120, 4310, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8250, 0, 3, 4226, 8130, 4328, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8280, 0, 3, 4232, 8140, 4346, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8310, 0, 3, 4238, 8150, 4364, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8340, 0, 3, 4244, 8160, 4382, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8370, 0, 3, 4256, 8170, 4418, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8400, 0, 3, 4262, 8180, 4436, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8430, 0, 3, 4268, 8190, 4454, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8460, 0, 3, 4274, 8200, 4472, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8490, 0, 3, 4280, 8210, 4490, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 8520, 0, 3, 4292, 8220, 1316, 1334,
                                                 4544, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8580, 0, 3, 4310, 8250, 1334, 1352,
                                                 4580, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8640, 0, 3, 4328, 8280, 1352, 1370,
                                                 4616, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8700, 0, 3, 4346, 8310, 1370, 1388,
                                                 4652, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8760, 0, 3, 4364, 8340, 1388, 1406,
                                                 4688, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8820, 0, 3, 4400, 8370, 1442, 1460,
                                                 4760, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8880, 0, 3, 4418, 8400, 1460, 1478,
                                                 4796, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8940, 0, 3, 4436, 8430, 1478, 1496,
                                                 4832, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9000, 0, 3, 4454, 8460, 1496, 1514,
                                                 4868, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9060, 0, 3, 4472, 8490, 1514, 1532,
                                                 4904, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9120, 0, 3, 4508, 8520, 1568, 1598,
                                                 4940, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9220, 0, 3, 4544, 8580, 1598, 1628,
                                                 5000, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9320, 0, 3, 4580, 8640, 1628, 1658,
                                                 5060, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9420, 0, 3, 4616, 8700, 1658, 1688,
                                                 5120, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9520, 0, 3, 4652, 8760, 1688, 1718,
                                                 5180, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9620, 0, 3, 4724, 8820, 1778, 1808,
                                                 5240, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9720, 0, 3, 4760, 8880, 1808, 1838,
                                                 5300, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9820, 0, 3, 4796, 8940, 1838, 1868,
                                                 5360, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9920, 0, 3, 4832, 9000, 1868, 1898,
                                                 5420, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10020, 0, 3, 4868, 9060, 1898, 1928,
                                                 5480, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10120, 0, 3, 8520, 8580, 5000, 9320,
                                                 1988, 2033, 5630, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10270, 0, 3, 8580, 8640, 5060, 9420,
                                                 2033, 2078, 5720, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10420, 0, 3, 8640, 8700, 5120, 9520,
                                                 2078, 2123, 5810, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10570, 0, 3, 8820, 8880, 5300, 9820,
                                                 2213, 2258, 5990, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10720, 0, 3, 8880, 8940, 5360, 9920,
                                                 2258, 2303, 6080, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10870, 0, 3, 8940, 9000, 5420, 10020,
                                                 2303, 2348, 6170, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11020, 0, 3, 9120, 9220, 5540, 10120,
                                                 2438, 2501, 6260, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11230, 0, 3, 9220, 9320, 5630, 10270,
                                                 2501, 2564, 6386, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11440, 0, 3, 9320, 9420, 5720, 10420,
                                                 2564, 2627, 6512, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11650, 0, 3, 9620, 9720, 5900, 10570,
                                                 2753, 2816, 6638, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11860, 0, 3, 9720, 9820, 5990, 10720,
                                                 2816, 2879, 6764, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12070, 0, 3, 9820, 9920, 6080, 10870,
                                                 2879, 2942, 6890, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12280, 0, 3, 10120, 10270, 6386, 11440,
                                                 3068, 3152, 7184, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12560, 0, 3, 10570, 10720, 6764, 12070,
                                                 3320, 3404, 7520, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12840, 0, 3, 11020, 11230, 7016, 12280,
                                                 3572, 3680, 7688, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 13200, 0, 3, 11650, 11860, 7352, 12560,
                                                 3896, 4004, 7904, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 13560, 12840, 720, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 14280, 13920, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 14280, 7, nmax);

    simdtrf::transform_f_inner(buffer, 14280, 13560, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 105 * nvalues, nvalues, buffer, 14280, 7, nmax);
}

}  // namespace simdt2ceri
