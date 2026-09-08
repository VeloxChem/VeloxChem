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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformLH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lh_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
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

    auto buffer = CSimdMatrix(8242, nvalues);

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

            compute_prim_hs_electron_repulsion_14(buffer, 230, 0, 89, 95, 158, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 242, 0, 95, 101, 167, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 254, 0, 101, 107, 176, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 266, 0, 107, 113, 185, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 278, 0, 113, 119, 194, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 290, 0, 119, 125, 203, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 302, 0, 125, 131, 212, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 314, 0, 131, 137, 221, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 326, 0, 149, 158, 242, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 341, 0, 158, 167, 254, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 356, 0, 167, 176, 266, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 371, 0, 176, 185, 278, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 386, 0, 185, 194, 290, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 401, 0, 194, 203, 302, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 416, 0, 203, 212, 314, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 431, 0, 230, 242, 341, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 446, 0, 242, 254, 356, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 464, 0, 254, 266, 371, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 482, 0, 266, 278, 386, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 500, 0, 278, 290, 401, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 518, 0, 290, 302, 416, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 536, 0, 326, 341, 446, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 554, 0, 341, 356, 464, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 572, 0, 356, 371, 482, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 590, 0, 371, 386, 500, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 608, 0, 386, 401, 518, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 626, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 629, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 632, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 635, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 638, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 641, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 644, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 647, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 650, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 653, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 656, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 659, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 662, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 665, 3, 12, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 668, 3, 13, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 671, 3, 14, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 674, 3, 15, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 677, 3, 16, 47, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 680, 3, 17, 50, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 683, 3, 20, 56, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 686, 3, 23, 59, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 689, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 692, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 695, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 698, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 701, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 704, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 707, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 710, 3, 47, 83, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 713, 3, 50, 86, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 716, 3, 59, 95, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 719, 3, 62, 101, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 722, 3, 65, 107, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 725, 3, 68, 113, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 728, 3, 71, 119, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 731, 3, 74, 125, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 734, 3, 77, 131, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 737, 3, 80, 137, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 740, 3, 83, 143, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 743, 3, 89, 149, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 746, 3, 95, 158, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 749, 3, 101, 167, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 752, 3, 107, 176, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 755, 3, 113, 185, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 764, 3, 119, 194, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 773, 3, 125, 203, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 782, 3, 131, 212, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 785, 3, 137, 221, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 788, 3, 158, 242, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 791, 3, 167, 254, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 794, 3, 176, 266, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 809, 3, 185, 278, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 824, 3, 194, 290, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 839, 3, 203, 302, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 848, 3, 212, 314, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 851, 3, 230, 326, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 854, 3, 242, 341, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 857, 3, 254, 356, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 875, 3, 266, 371, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 893, 3, 278, 386, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 911, 3, 290, 401, ncols, p);

            compute_prim_ip_electron_repulsion_17(buffer, 929, 3, 302, 416, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 938, 3, 341, 446, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 959, 3, 356, 464, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 980, 3, 371, 482, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1001, 3, 386, 500, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1022, 3, 401, 518, ncols, p);

            compute_prim_lp_electron_repulsion_2(buffer, 1043, 3, 431, 536, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1067, 3, 446, 554, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1091, 3, 464, 572, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1115, 3, 482, 590, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1139, 3, 500, 608, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1163, 3, 8, 9, 629, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1166, 3, 9, 10, 632, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1169, 3, 10, 11, 635, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1172, 3, 11, 12, 638, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1175, 3, 12, 13, 641, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1178, 3, 13, 14, 644, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1181, 3, 14, 15, 647, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1184, 3, 15, 16, 650, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1187, 3, 16, 17, 653, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1190, 0, 626, 1163, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1193, 0, 629, 1166, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1196, 0, 632, 1169, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1199, 0, 635, 1172, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1202, 0, 638, 1175, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1205, 0, 641, 1178, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1208, 0, 644, 1181, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1211, 0, 647, 1184, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1214, 0, 650, 1187, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1217, 3, 656, 56, 59, 689, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1220, 3, 659, 59, 62, 692, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1223, 3, 662, 62, 65, 695, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1226, 3, 665, 65, 68, 698, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1229, 3, 668, 68, 71, 701, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1232, 3, 671, 71, 74, 704, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1235, 3, 674, 74, 77, 707, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1238, 3, 677, 77, 80, 710, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1241, 3, 680, 80, 83, 713, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1244, 0, 3, 689, 1220, 89, 95, 719, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1253, 0, 3, 692, 1223, 95, 101, 722, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1262, 0, 3, 695, 1226, 101, 107, 725, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1271, 0, 3, 698, 1229, 107, 113, 728, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1280, 0, 3, 701, 1232, 113, 119, 731, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1289, 0, 3, 704, 1235, 119, 125, 734, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1298, 0, 3, 707, 1238, 125, 131, 737, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1307, 0, 3, 710, 1241, 131, 137, 740, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1316, 0, 3, 1217, 1220, 719, 1253, 149, 158, 749, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1331, 0, 3, 1220, 1223, 722, 1262, 158, 167, 752, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_39(buffer, 1346, 0, 3, 1223, 1226, 725, 1271, 167, 176, 755, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1361, 0, 3, 1226, 1229, 728, 1280, 176, 185, 764, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1382, 0, 3, 1229, 1232, 731, 1289, 185, 194, 773, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1403, 0, 3, 1232, 1235, 734, 1298, 194, 203, 782, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1418, 0, 3, 1235, 1238, 737, 1307, 203, 212, 785, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1433, 0, 3, 1244, 1253, 749, 1331, 230, 242, 791, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1454, 0, 3, 1253, 1262, 752, 1346, 242, 254, 794, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1475, 0, 3, 1262, 1271, 755, 1361, 254, 266, 809, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1508, 0, 3, 1271, 1280, 764, 1382, 266, 278, 824, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_45(buffer, 1541, 0, 3, 1280, 1289, 773, 1403, 278, 290, 839, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1568, 0, 3, 1289, 1298, 782, 1418, 290, 302, 848, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_33(buffer, 1589, 0, 3, 1316, 1331, 791, 1454, 326, 341, 857, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 1616, 0, 3, 1331, 1346, 794, 1475, 341, 356, 875, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_34(buffer, 1658, 0, 3, 1346, 1361, 809, 1508, 356, 371, 893, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_35(buffer, 1703, 0, 3, 1361, 1382, 824, 1541, 371, 386, 911, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_36(buffer, 1745, 0, 3, 1382, 1403, 839, 1568, 386, 401, 929, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_15(buffer, 1778, 0, 3, 1433, 1454, 857, 1616, 431, 446, 959, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_16(buffer, 1826, 0, 3, 1454, 1475, 875, 1658, 446, 464, 980, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_17(buffer, 1877, 0, 3, 1475, 1508, 893, 1703, 464, 482, 1001, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_18(buffer, 1931, 0, 3, 1508, 1541, 911, 1745, 482, 500, 1022, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_5(buffer, 1982, 0, 3, 1589, 1616, 959, 1826, 536, 554, 1091, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_6(buffer, 2036, 0, 3, 1616, 1658, 980, 1877, 554, 572, 1115, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_7(buffer, 2090, 0, 3, 1658, 1703, 1001, 1931, 572, 590, 1139, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2144, 3, 626, 629, 1166, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2147, 3, 629, 632, 1169, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2150, 3, 632, 635, 1172, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2153, 3, 635, 638, 1175, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2156, 3, 638, 641, 1178, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2159, 3, 641, 644, 1181, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2162, 3, 644, 647, 1184, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2165, 3, 647, 650, 1187, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2168, 0, 1166, 2147, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2171, 0, 1169, 2150, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2174, 0, 1172, 2153, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2177, 0, 1175, 2156, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2180, 0, 1178, 2159, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2183, 0, 1181, 2162, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2186, 0, 1184, 2165, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2189, 3, 1190, 683, 686, 1217, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2192, 3, 1193, 686, 689, 1220, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2195, 3, 1196, 689, 692, 1223, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2201, 3, 1199, 692, 695, 1226, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2207, 3, 1202, 695, 698, 1229, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2213, 3, 1205, 698, 701, 1232, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2219, 3, 1208, 701, 704, 1235, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2225, 3, 1211, 704, 707, 1238, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2231, 3, 1214, 707, 710, 1241, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_31(buffer, 2234, 0, 3, 1220, 2195, 716, 719, 1253, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 2243, 0, 3, 1223, 2201, 719, 722, 1262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2255, 0, 3, 1226, 2207, 722, 725, 1271, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2270, 0, 3, 1229, 2213, 725, 728, 1280, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2285, 0, 3, 1232, 2219, 728, 731, 1289, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2300, 0, 3, 1235, 2225, 731, 734, 1298, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2315, 0, 3, 1238, 2231, 734, 737, 1307, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_29(buffer, 2324, 0, 3, 2189, 2192, 1244, 2234, 743, 746, 1316, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_57(buffer, 2342, 0, 3, 2192, 2195, 1253, 2243, 746, 749, 1331, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_58(buffer, 2360, 0, 3, 2195, 2201, 1262, 2255, 749, 752, 1346, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_59(buffer, 2384, 0, 3, 2201, 2207, 1271, 2270, 752, 755, 1361, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 2408, 0, 3, 2207, 2213, 1280, 2285, 755, 764, 1382, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_61(buffer, 2438, 0, 3, 2213, 2219, 1289, 2300, 764, 773, 1403, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_62(buffer, 2462, 0, 3, 2219, 2225, 1298, 2315, 773, 782, 1418, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_47(buffer, 2480, 0, 3, 2234, 2243, 1331, 2360, 788, 791, 1454, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_48(buffer, 2507, 0, 3, 2243, 2255, 1346, 2384, 791, 794, 1475, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_49(buffer, 2540, 0, 3, 2255, 2270, 1361, 2408, 794, 809, 1508, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_50(buffer, 2600, 0, 3, 2270, 2285, 1382, 2438, 809, 824, 1541, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_51(buffer, 2648, 0, 3, 2285, 2300, 1403, 2462, 824, 839, 1568, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_29(buffer, 2678, 0, 3, 2324, 2342, 1433, 2480, 851, 854, 1589, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_30(buffer, 2714, 0, 3, 2342, 2360, 1454, 2507, 854, 857, 1616, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_31(buffer, 2750, 0, 3, 2360, 2384, 1475, 2540, 857, 875, 1658, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_32(buffer, 2870, 0, 3, 2384, 2408, 1508, 2600, 875, 893, 1703, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_33(buffer, 2969, 0, 3, 2408, 2438, 1541, 2648, 893, 911, 1745, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_13(buffer, 3029, 0, 3, 2480, 2507, 1616, 2750, 938, 959, 1826, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_14(buffer, 3125, 0, 3, 2507, 2540, 1658, 2870, 959, 980, 1877, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_15(buffer, 3269, 0, 3, 2540, 2600, 1703, 2969, 980, 1001, 1931, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_2(buffer, 3389, 0, 3, 2678, 2714, 1778, 3029, 1043, 1067, 1982, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_3(buffer, 3497, 0, 3, 2714, 2750, 1826, 3125, 1067, 1091, 2036, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_4(buffer, 3605, 0, 3, 2750, 2870, 1877, 3269, 1091, 1115, 2090, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3773, 3, 1163, 1166, 2147, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3776, 3, 1166, 1169, 2150, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3779, 3, 1169, 1172, 2153, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3782, 3, 1172, 1175, 2156, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3785, 3, 1175, 1178, 2159, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3788, 3, 1178, 1181, 2162, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3791, 3, 1181, 1184, 2165, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3794, 0, 2144, 3773, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3797, 0, 2147, 3776, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3800, 0, 2150, 3779, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3803, 0, 2153, 3782, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3806, 0, 2156, 3785, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3809, 0, 2159, 3788, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3812, 0, 2162, 3791, ncols, p);

            compute_prim_dg_electron_repulsion_27(buffer, 3815, 3, 2168, 1217, 1220, 2195, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 3818, 3, 2171, 1220, 1223, 2201, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 3824, 3, 2174, 1223, 1226, 2207, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 3830, 3, 2177, 1226, 1229, 2213, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 3836, 3, 2180, 1229, 1232, 2219, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 3842, 3, 2183, 1232, 1235, 2225, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3848, 3, 2186, 1235, 1238, 2231, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_30(buffer, 3851, 0, 3, 2195, 3818, 1244, 1253, 2243, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_53(buffer, 3866, 0, 3, 2201, 3824, 1253, 1262, 2255, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_45(buffer, 3881, 0, 3, 2207, 3830, 1262, 1271, 2270, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_45(buffer, 3899, 0, 3, 2213, 3836, 1271, 1280, 2285, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_45(buffer, 3917, 0, 3, 2219, 3842, 1280, 1289, 2300, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_54(buffer, 3935, 0, 3, 2225, 3848, 1289, 1298, 2315, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_57(buffer, 3944, 0, 3, 3815, 3818, 2243, 3866, 1316, 1331, 2360, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_58(buffer, 3965, 0, 3, 3818, 3824, 2255, 3881, 1331, 1346, 2384, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_59(buffer, 3995, 0, 3, 3824, 3830, 2270, 3899, 1346, 1361, 2408, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_60(buffer, 4025, 0, 3, 3830, 3836, 2285, 3917, 1361, 1382, 2438, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_61(buffer, 4055, 0, 3, 3836, 3842, 2300, 3935, 1382, 1403, 2462, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_44(buffer, 4076, 0, 3, 3851, 3866, 2360, 3965, 1433, 1454, 2507, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_45(buffer, 4121, 0, 3, 3866, 3881, 2384, 3995, 1454, 1475, 2540, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_46(buffer, 4166, 0, 3, 3881, 3899, 2408, 4025, 1475, 1508, 2600, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_47(buffer, 4226, 0, 3, 3899, 3917, 2438, 4055, 1508, 1541, 2648, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_26(buffer, 4265, 0, 3, 3944, 3965, 2507, 4121, 1589, 1616, 2750, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_27(buffer, 4328, 0, 3, 3965, 3995, 2540, 4166, 1616, 1658, 2870, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_28(buffer, 4474, 0, 3, 3995, 4025, 2600, 4226, 1658, 1703, 2969, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_11(buffer, 4561, 0, 3, 4076, 4121, 2750, 4328, 1778, 1826, 3125, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_12(buffer, 4907, 0, 3, 4121, 4166, 2870, 4474, 1826, 1877, 3269, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_1(buffer, 5116, 0, 3, 4265, 4328, 3125, 4907, 1982, 2036, 3605, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_9(buffer, 5554, 3, 3794, 2189, 2192, 3815, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_28(buffer, 5557, 3, 3797, 2192, 2195, 3818, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 5560, 3, 3800, 2195, 2201, 3824, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 5563, 3, 3803, 2201, 2207, 3830, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 5566, 3, 3806, 2207, 2213, 3836, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 5569, 3, 3809, 2213, 2219, 3842, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 5572, 3, 3812, 2219, 2225, 3848, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_27(buffer, 5575, 0, 3, 3818, 5560, 2234, 2243, 3866, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_55(buffer, 5584, 0, 3, 3824, 5563, 2243, 2255, 3881, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_43(buffer, 5593, 0, 3, 3830, 5566, 2255, 2270, 3899, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_43(buffer, 5602, 0, 3, 3836, 5569, 2270, 2285, 3917, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_56(buffer, 5611, 0, 3, 3842, 5572, 2285, 2300, 3935, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_18(buffer, 5620, 0, 3, 5554, 5557, 3851, 5575, 2324, 2342, 3944, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_48(buffer, 5641, 0, 3, 5557, 5560, 3866, 5584, 2342, 2360, 3965, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_49(buffer, 5662, 0, 3, 5560, 5563, 3881, 5593, 2360, 2384, 3995, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_50(buffer, 5683, 0, 3, 5563, 5566, 3899, 5602, 2384, 2408, 4025, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_51(buffer, 5704, 0, 3, 5566, 5569, 3917, 5611, 2408, 2438, 4055, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_33(buffer, 5725, 0, 3, 5575, 5584, 3965, 5662, 2480, 2507, 4121, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_34(buffer, 5767, 0, 3, 5584, 5593, 3995, 5683, 2507, 2540, 4166, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_35(buffer, 5809, 0, 3, 5593, 5602, 4025, 5704, 2540, 2600, 4226, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_18(buffer, 5851, 0, 3, 5620, 5641, 4076, 5725, 2678, 2714, 4265, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_19(buffer, 5923, 0, 3, 5641, 5662, 4121, 5767, 2714, 2750, 4328, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_20(buffer, 5995, 0, 3, 5662, 5683, 4166, 5809, 2750, 2870, 4474, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_8(buffer, 6091, 0, 3, 5725, 5767, 4328, 5995, 3029, 3125, 4907, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 6352, 0, 3, 5851, 5923, 4561, 6091, 3389, 3497, 5116, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 7297, 6352, 945, ncols);
        }
    }

    simdtrf::transform_lh(values, nvalues, buffer, 7297, nmax);
}

}  // namespace simdt2ceri
