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


#include "SimdElectronRepulsionRecLH.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
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
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 34618, 33178, 945, nvalues);

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
                                            10, 11, 12, 13}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 7, 8, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 8, 9, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 9, 10, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 10, 11, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 11, 12, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 12, 13, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 13, 14, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 14, 15, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 15, 16, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 16, 17, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 20, 23, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 23, 26, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 26, 29, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 29, 32, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 32, 35, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 35, 38, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 38, 41, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 41, 44, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 44, 47, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 47, 50, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 222, 0, 56, 62, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 237, 0, 62, 68, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 252, 0, 68, 74, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 267, 0, 74, 80, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 282, 0, 80, 86, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 297, 0, 86, 92, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 312, 0, 92, 98, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 327, 0, 98, 104, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 342, 0, 104, 110, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 357, 0, 122, 132, 237, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 378, 0, 132, 142, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 399, 0, 142, 152, 267, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 420, 0, 152, 162, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 441, 0, 162, 172, 297, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 462, 0, 172, 182, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 483, 0, 182, 192, 327, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 504, 0, 192, 202, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 222, 237, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 553, 0, 237, 252, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 252, 267, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 267, 282, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 282, 297, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 297, 312, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 312, 327, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 721, 0, 357, 378, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 757, 0, 378, 399, 581, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 793, 0, 399, 420, 609, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 829, 0, 420, 441, 637, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 865, 0, 441, 462, 665, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 901, 0, 462, 483, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 937, 0, 525, 553, 757, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 982, 0, 553, 581, 793, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1027, 0, 581, 609, 829, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1072, 0, 609, 637, 865, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1117, 0, 637, 665, 901, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1162, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1165, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1168, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1171, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1174, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1177, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1180, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1183, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1186, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1189, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1192, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1195, 3, 8, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1204, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1213, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1222, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1231, 3, 12, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1240, 3, 13, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1249, 3, 14, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1258, 3, 15, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1267, 3, 16, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1276, 3, 17, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1285, 3, 18, 53, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1294, 0, 3, 20, 1195, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1312, 0, 3, 23, 1204, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1330, 0, 3, 26, 1213, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1348, 0, 3, 29, 1222, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1366, 0, 3, 32, 1231, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1384, 0, 3, 35, 1240, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1402, 0, 3, 38, 1249, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1420, 0, 3, 41, 1258, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1438, 0, 3, 44, 1267, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1456, 0, 3, 47, 1276, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1474, 0, 3, 50, 1285, 116, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1492, 0, 3, 62, 1330, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1522, 0, 3, 68, 1348, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1552, 0, 3, 74, 1366, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1582, 0, 3, 80, 1384, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1612, 0, 3, 86, 1402, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1642, 0, 3, 92, 1420, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1672, 0, 3, 98, 1438, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1702, 0, 3, 104, 1456, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1732, 0, 3, 110, 1474, 212, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1762, 0, 3, 122, 1492, 222, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1807, 0, 3, 132, 1522, 237, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1852, 0, 3, 142, 1552, 252, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1897, 0, 3, 152, 1582, 267, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1942, 0, 3, 162, 1612, 282, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1987, 0, 3, 172, 1642, 297, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2032, 0, 3, 182, 1672, 312, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2077, 0, 3, 192, 1702, 327, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2122, 0, 3, 202, 1732, 342, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2167, 0, 3, 237, 1852, 378, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2230, 0, 3, 252, 1897, 399, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2293, 0, 3, 267, 1942, 420, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2356, 0, 3, 282, 1987, 441, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2419, 0, 3, 297, 2032, 462, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2482, 0, 3, 312, 2077, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2545, 0, 3, 327, 2122, 504, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2608, 0, 3, 357, 2167, 525, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2692, 0, 3, 378, 2230, 553, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2776, 0, 3, 399, 2293, 581, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2860, 0, 3, 420, 2356, 609, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2944, 0, 3, 441, 2419, 637, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3028, 0, 3, 462, 2482, 665, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3112, 0, 3, 483, 2545, 693, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3196, 0, 3, 553, 2776, 757, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3304, 0, 3, 581, 2860, 793, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3412, 0, 3, 609, 2944, 829, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3520, 0, 3, 637, 3028, 865, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3628, 0, 3, 665, 3112, 901, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3736, 0, 3, 721, 3196, 937, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3871, 0, 3, 757, 3304, 982, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4006, 0, 3, 793, 3412, 1027, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4141, 0, 3, 829, 3520, 1072, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4276, 0, 3, 865, 3628, 1117, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4411, 3, 8, 9, 1165, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4417, 3, 9, 10, 1168, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4423, 3, 10, 11, 1171, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4429, 3, 11, 12, 1174, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4435, 3, 12, 13, 1177, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4441, 3, 13, 14, 1180, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4447, 3, 14, 15, 1183, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4453, 3, 15, 16, 1186, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4459, 3, 16, 17, 1189, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4465, 3, 17, 18, 1192, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4471, 0, 3, 1162, 4411, 1204, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4489, 0, 3, 1165, 4417, 1213, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4507, 0, 3, 1168, 4423, 1222, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4525, 0, 3, 1171, 4429, 1231, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4543, 0, 3, 1174, 4435, 1240, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4561, 0, 3, 1177, 4441, 1249, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4579, 0, 3, 1180, 4447, 1258, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4597, 0, 3, 1183, 4453, 1267, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4615, 0, 3, 1186, 4459, 1276, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4633, 0, 3, 1189, 4465, 1285, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4651, 0, 3, 1204, 4489, 56, 62, 1330,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4687, 0, 3, 1213, 4507, 62, 68, 1348,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4723, 0, 3, 1222, 4525, 68, 74, 1366,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4759, 0, 3, 1231, 4543, 74, 80, 1384,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4795, 0, 3, 1240, 4561, 80, 86, 1402,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4831, 0, 3, 1249, 4579, 86, 92, 1420,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4867, 0, 3, 1258, 4597, 92, 98, 1438,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4903, 0, 3, 1267, 4615, 98, 104, 1456,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4939, 0, 3, 1276, 4633, 104, 110, 1474,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4975, 0, 3, 1330, 4687, 122, 132, 1522,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5035, 0, 3, 1348, 4723, 132, 142, 1552,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5095, 0, 3, 1366, 4759, 142, 152, 1582,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5155, 0, 3, 1384, 4795, 152, 162, 1612,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5215, 0, 3, 1402, 4831, 162, 172, 1642,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5275, 0, 3, 1420, 4867, 172, 182, 1672,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5335, 0, 3, 1438, 4903, 182, 192, 1702,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5395, 0, 3, 1456, 4939, 192, 202, 1732,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5455, 0, 3, 4651, 4687, 1522, 5035, 222,
                                                 237, 1852, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5545, 0, 3, 4687, 4723, 1552, 5095, 237,
                                                 252, 1897, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5635, 0, 3, 4723, 4759, 1582, 5155, 252,
                                                 267, 1942, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5725, 0, 3, 4759, 4795, 1612, 5215, 267,
                                                 282, 1987, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5815, 0, 3, 4795, 4831, 1642, 5275, 282,
                                                 297, 2032, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5905, 0, 3, 4831, 4867, 1672, 5335, 297,
                                                 312, 2077, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5995, 0, 3, 4867, 4903, 1702, 5395, 312,
                                                 327, 2122, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6085, 0, 3, 4975, 5035, 1852, 5545, 357,
                                                 378, 2230, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6211, 0, 3, 5035, 5095, 1897, 5635, 378,
                                                 399, 2293, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6337, 0, 3, 5095, 5155, 1942, 5725, 399,
                                                 420, 2356, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6463, 0, 3, 5155, 5215, 1987, 5815, 420,
                                                 441, 2419, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6589, 0, 3, 5215, 5275, 2032, 5905, 441,
                                                 462, 2482, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6715, 0, 3, 5275, 5335, 2077, 5995, 462,
                                                 483, 2545, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6841, 0, 3, 5455, 5545, 2230, 6211, 525,
                                                 553, 2776, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7009, 0, 3, 5545, 5635, 2293, 6337, 553,
                                                 581, 2860, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7177, 0, 3, 5635, 5725, 2356, 6463, 581,
                                                 609, 2944, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7345, 0, 3, 5725, 5815, 2419, 6589, 609,
                                                 637, 3028, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7513, 0, 3, 5815, 5905, 2482, 6715, 637,
                                                 665, 3112, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7681, 0, 3, 6085, 6211, 2776, 7009, 721,
                                                 757, 3304, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7897, 0, 3, 6211, 6337, 2860, 7177, 757,
                                                 793, 3412, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8113, 0, 3, 6337, 6463, 2944, 7345, 793,
                                                 829, 3520, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8329, 0, 3, 6463, 6589, 3028, 7513, 829,
                                                 865, 3628, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 8545, 0, 3, 6841, 7009, 3304, 7897, 937,
                                                 982, 4006, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 8815, 0, 3, 7009, 7177, 3412, 8113, 982,
                                                 1027, 4141, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 9085, 0, 3, 7177, 7345, 3520, 8329,
                                                 1027, 1072, 4276, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9355, 3, 1162, 1165, 4417, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9365, 3, 1165, 1168, 4423, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9375, 3, 1168, 1171, 4429, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9385, 3, 1171, 1174, 4435, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9395, 3, 1174, 1177, 4441, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9405, 3, 1177, 1180, 4447, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9415, 3, 1180, 1183, 4453, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9425, 3, 1183, 1186, 4459, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9435, 3, 1186, 1189, 4465, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 9445, 0, 3, 4411, 9355, 4489, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9475, 0, 3, 4417, 9365, 4507, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9505, 0, 3, 4423, 9375, 4525, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9535, 0, 3, 4429, 9385, 4543, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9565, 0, 3, 4435, 9395, 4561, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9595, 0, 3, 4441, 9405, 4579, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9625, 0, 3, 4447, 9415, 4597, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9655, 0, 3, 4453, 9425, 4615, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9685, 0, 3, 4459, 9435, 4633, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 9715, 0, 3, 4471, 9445, 1294, 1312,
                                                 4651, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9775, 0, 3, 4489, 9475, 1312, 1330,
                                                 4687, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9835, 0, 3, 4507, 9505, 1330, 1348,
                                                 4723, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9895, 0, 3, 4525, 9535, 1348, 1366,
                                                 4759, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9955, 0, 3, 4543, 9565, 1366, 1384,
                                                 4795, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10015, 0, 3, 4561, 9595, 1384, 1402,
                                                 4831, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10075, 0, 3, 4579, 9625, 1402, 1420,
                                                 4867, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10135, 0, 3, 4597, 9655, 1420, 1438,
                                                 4903, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10195, 0, 3, 4615, 9685, 1438, 1456,
                                                 4939, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10255, 0, 3, 4687, 9835, 1492, 1522,
                                                 5035, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10355, 0, 3, 4723, 9895, 1522, 1552,
                                                 5095, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10455, 0, 3, 4759, 9955, 1552, 1582,
                                                 5155, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10555, 0, 3, 4795, 10015, 1582, 1612,
                                                 5215, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10655, 0, 3, 4831, 10075, 1612, 1642,
                                                 5275, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10755, 0, 3, 4867, 10135, 1642, 1672,
                                                 5335, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10855, 0, 3, 4903, 10195, 1672, 1702,
                                                 5395, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10955, 0, 3, 9715, 9775, 4975, 10255,
                                                 1762, 1807, 5455, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11105, 0, 3, 9775, 9835, 5035, 10355,
                                                 1807, 1852, 5545, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11255, 0, 3, 9835, 9895, 5095, 10455,
                                                 1852, 1897, 5635, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11405, 0, 3, 9895, 9955, 5155, 10555,
                                                 1897, 1942, 5725, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11555, 0, 3, 9955, 10015, 5215, 10655,
                                                 1942, 1987, 5815, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11705, 0, 3, 10015, 10075, 5275, 10755,
                                                 1987, 2032, 5905, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11855, 0, 3, 10075, 10135, 5335, 10855,
                                                 2032, 2077, 5995, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12005, 0, 3, 10255, 10355, 5545, 11255,
                                                 2167, 2230, 6211, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12215, 0, 3, 10355, 10455, 5635, 11405,
                                                 2230, 2293, 6337, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12425, 0, 3, 10455, 10555, 5725, 11555,
                                                 2293, 2356, 6463, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12635, 0, 3, 10555, 10655, 5815, 11705,
                                                 2356, 2419, 6589, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12845, 0, 3, 10655, 10755, 5905, 11855,
                                                 2419, 2482, 6715, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13055, 0, 3, 10955, 11105, 6085, 12005,
                                                 2608, 2692, 6841, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13335, 0, 3, 11105, 11255, 6211, 12215,
                                                 2692, 2776, 7009, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13615, 0, 3, 11255, 11405, 6337, 12425,
                                                 2776, 2860, 7177, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13895, 0, 3, 11405, 11555, 6463, 12635,
                                                 2860, 2944, 7345, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14175, 0, 3, 11555, 11705, 6589, 12845,
                                                 2944, 3028, 7513, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 14455, 0, 3, 12005, 12215, 7009, 13615,
                                                 3196, 3304, 7897, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 14815, 0, 3, 12215, 12425, 7177, 13895,
                                                 3304, 3412, 8113, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 15175, 0, 3, 12425, 12635, 7345, 14175,
                                                 3412, 3520, 8329, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 15535, 0, 3, 13055, 13335, 7681, 14455,
                                                 3736, 3871, 8545, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 15985, 0, 3, 13335, 13615, 7897, 14815,
                                                 3871, 4006, 8815, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 16435, 0, 3, 13615, 13895, 8113, 15175,
                                                 4006, 4141, 9085, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16885, 3, 4411, 4417, 9365, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16900, 3, 4417, 4423, 9375, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16915, 3, 4423, 4429, 9385, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16930, 3, 4429, 4435, 9395, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16945, 3, 4435, 4441, 9405, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16960, 3, 4441, 4447, 9415, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16975, 3, 4447, 4453, 9425, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16990, 3, 4453, 4459, 9435, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 17005, 0, 3, 9355, 16885, 9475, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17050, 0, 3, 9365, 16900, 9505, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17095, 0, 3, 9375, 16915, 9535, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17140, 0, 3, 9385, 16930, 9565, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17185, 0, 3, 9395, 16945, 9595, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17230, 0, 3, 9405, 16960, 9625, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17275, 0, 3, 9415, 16975, 9655, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17320, 0, 3, 9425, 16990, 9685, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 17365, 0, 3, 9475, 17050, 4651, 4687,
                                                 9835, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17455, 0, 3, 9505, 17095, 4687, 4723,
                                                 9895, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17545, 0, 3, 9535, 17140, 4723, 4759,
                                                 9955, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17635, 0, 3, 9565, 17185, 4759, 4795,
                                                 10015, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17725, 0, 3, 9595, 17230, 4795, 4831,
                                                 10075, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17815, 0, 3, 9625, 17275, 4831, 4867,
                                                 10135, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17905, 0, 3, 9655, 17320, 4867, 4903,
                                                 10195, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17995, 0, 3, 9835, 17455, 4975, 5035,
                                                 10355, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18145, 0, 3, 9895, 17545, 5035, 5095,
                                                 10455, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18295, 0, 3, 9955, 17635, 5095, 5155,
                                                 10555, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18445, 0, 3, 10015, 17725, 5155, 5215,
                                                 10655, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18595, 0, 3, 10075, 17815, 5215, 5275,
                                                 10755, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18745, 0, 3, 10135, 17905, 5275, 5335,
                                                 10855, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18895, 0, 3, 17365, 17455, 10355, 18145,
                                                 5455, 5545, 11255, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19120, 0, 3, 17455, 17545, 10455, 18295,
                                                 5545, 5635, 11405, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19345, 0, 3, 17545, 17635, 10555, 18445,
                                                 5635, 5725, 11555, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19570, 0, 3, 17635, 17725, 10655, 18595,
                                                 5725, 5815, 11705, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19795, 0, 3, 17725, 17815, 10755, 18745,
                                                 5815, 5905, 11855, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20020, 0, 3, 17995, 18145, 11255, 19120,
                                                 6085, 6211, 12215, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20335, 0, 3, 18145, 18295, 11405, 19345,
                                                 6211, 6337, 12425, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20650, 0, 3, 18295, 18445, 11555, 19570,
                                                 6337, 6463, 12635, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20965, 0, 3, 18445, 18595, 11705, 19795,
                                                 6463, 6589, 12845, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 21280, 0, 3, 18895, 19120, 12215, 20335,
                                                 6841, 7009, 13615, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 21700, 0, 3, 19120, 19345, 12425, 20650,
                                                 7009, 7177, 13895, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 22120, 0, 3, 19345, 19570, 12635, 20965,
                                                 7177, 7345, 14175, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 22540, 0, 3, 20020, 20335, 13615, 21700,
                                                 7681, 7897, 14815, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 23080, 0, 3, 20335, 20650, 13895, 22120,
                                                 7897, 8113, 15175, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 23620, 0, 3, 21280, 21700, 14815, 23080,
                                                 8545, 8815, 16435, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24295, 3, 9355, 9365, 16900, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24316, 3, 9365, 9375, 16915, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24337, 3, 9375, 9385, 16930, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24358, 3, 9385, 9395, 16945, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24379, 3, 9395, 9405, 16960, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24400, 3, 9405, 9415, 16975, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24421, 3, 9415, 9425, 16990, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 24442, 0, 3, 16885, 24295, 17050, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24505, 0, 3, 16900, 24316, 17095, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24568, 0, 3, 16915, 24337, 17140, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24631, 0, 3, 16930, 24358, 17185, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24694, 0, 3, 16945, 24379, 17230, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24757, 0, 3, 16960, 24400, 17275, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 24820, 0, 3, 16975, 24421, 17320, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 24883, 0, 3, 17005, 24442, 9715, 9775,
                                                 17365, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25009, 0, 3, 17050, 24505, 9775, 9835,
                                                 17455, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25135, 0, 3, 17095, 24568, 9835, 9895,
                                                 17545, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25261, 0, 3, 17140, 24631, 9895, 9955,
                                                 17635, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25387, 0, 3, 17185, 24694, 9955, 10015,
                                                 17725, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25513, 0, 3, 17230, 24757, 10015, 10075,
                                                 17815, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25639, 0, 3, 17275, 24820, 10075, 10135,
                                                 17905, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 25765, 0, 3, 17455, 25135, 10255, 10355,
                                                 18145, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 25975, 0, 3, 17545, 25261, 10355, 10455,
                                                 18295, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26185, 0, 3, 17635, 25387, 10455, 10555,
                                                 18445, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26395, 0, 3, 17725, 25513, 10555, 10655,
                                                 18595, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26605, 0, 3, 17815, 25639, 10655, 10755,
                                                 18745, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26815, 0, 3, 24883, 25009, 17995, 25765,
                                                 10955, 11105, 18895, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 27130, 0, 3, 25009, 25135, 18145, 25975,
                                                 11105, 11255, 19120, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 27445, 0, 3, 25135, 25261, 18295, 26185,
                                                 11255, 11405, 19345, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 27760, 0, 3, 25261, 25387, 18445, 26395,
                                                 11405, 11555, 19570, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 28075, 0, 3, 25387, 25513, 18595, 26605,
                                                 11555, 11705, 19795, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 28390, 0, 3, 25765, 25975, 19120, 27445,
                                                 12005, 12215, 20335, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 28831, 0, 3, 25975, 26185, 19345, 27760,
                                                 12215, 12425, 20650, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 29272, 0, 3, 26185, 26395, 19570, 28075,
                                                 12425, 12635, 20965, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 29713, 0, 3, 26815, 27130, 20020, 28390,
                                                 13055, 13335, 21280, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 30301, 0, 3, 27130, 27445, 20335, 28831,
                                                 13335, 13615, 21700, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 30889, 0, 3, 27445, 27760, 20650, 29272,
                                                 13615, 13895, 22120, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 31477, 0, 3, 28390, 28831, 21700, 30889,
                                                 14455, 14815, 23080, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 32233, 0, 3, 29713, 30301, 22540, 31477,
                                                 15535, 15985, 23620, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 33178, 32233, 945, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 34123, 33178, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 34123, 11, nmax);
}

}  // namespace simdt2ceri
