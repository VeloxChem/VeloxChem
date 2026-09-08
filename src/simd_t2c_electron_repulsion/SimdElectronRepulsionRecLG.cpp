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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformLG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lg_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
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

    auto buffer = CSimdMatrix(5263, nvalues);

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

            compute_prim_ds_electron_repulsion_1(buffer, 53, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 68, 0, 12, 13, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 71, 0, 13, 14, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 74, 0, 14, 15, 41, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 77, 0, 15, 16, 44, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 80, 0, 16, 17, 47, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 83, 0, 17, 18, 50, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 86, 0, 20, 23, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 92, 0, 23, 26, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 98, 0, 26, 29, 65, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 104, 0, 29, 32, 68, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 110, 0, 32, 35, 71, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 116, 0, 35, 38, 74, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 122, 0, 38, 41, 77, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 128, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 134, 0, 44, 47, 83, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 140, 0, 53, 56, 86, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 149, 0, 56, 59, 92, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 158, 0, 59, 62, 98, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 167, 0, 62, 65, 104, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 176, 0, 65, 68, 110, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 185, 0, 68, 71, 116, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 194, 0, 71, 74, 122, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 203, 0, 74, 77, 128, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 212, 0, 77, 80, 134, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 221, 0, 86, 92, 158, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 233, 0, 92, 98, 167, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 245, 0, 98, 104, 176, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 257, 0, 104, 110, 185, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 269, 0, 110, 116, 194, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 281, 0, 116, 122, 203, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 293, 0, 122, 128, 212, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 305, 0, 140, 149, 221, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 320, 0, 149, 158, 233, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 335, 0, 158, 167, 245, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 350, 0, 167, 176, 257, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 365, 0, 176, 185, 269, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 380, 0, 185, 194, 281, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 395, 0, 194, 203, 293, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 410, 0, 221, 233, 335, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 428, 0, 233, 245, 350, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 446, 0, 245, 257, 365, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 464, 0, 257, 269, 380, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 482, 0, 269, 281, 395, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 500, 0, 305, 320, 410, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 518, 0, 320, 335, 428, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 536, 0, 335, 350, 446, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 554, 0, 350, 365, 464, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 572, 0, 365, 380, 482, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 590, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 593, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 596, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 599, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 602, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 605, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 608, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 611, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 614, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 617, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 620, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 623, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 626, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 629, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 632, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 635, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 638, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 641, 3, 17, 47, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 644, 3, 23, 59, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 647, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 650, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 653, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 656, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 659, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 662, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 665, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 668, 3, 47, 83, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 671, 3, 59, 92, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 674, 3, 62, 98, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 677, 3, 65, 104, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 680, 3, 68, 110, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 683, 3, 71, 116, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 686, 3, 74, 122, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 689, 3, 77, 128, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 692, 3, 80, 134, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 695, 3, 92, 158, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 698, 3, 98, 167, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 701, 3, 104, 176, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 704, 3, 110, 185, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 713, 3, 116, 194, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 722, 3, 122, 203, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 725, 3, 128, 212, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 728, 3, 158, 233, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 731, 3, 167, 245, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 734, 3, 176, 257, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 749, 3, 185, 269, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 764, 3, 194, 281, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 773, 3, 203, 293, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 776, 3, 233, 335, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 779, 3, 245, 350, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 797, 3, 257, 365, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 815, 3, 269, 380, ncols, p);

            compute_prim_ip_electron_repulsion_17(buffer, 833, 3, 281, 395, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 842, 3, 335, 428, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 863, 3, 350, 446, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 884, 3, 365, 464, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 905, 3, 380, 482, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 926, 3, 428, 536, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 950, 3, 446, 554, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 974, 3, 464, 572, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 998, 3, 9, 10, 593, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1001, 3, 10, 11, 596, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1004, 3, 11, 12, 599, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1007, 3, 12, 13, 602, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1010, 3, 13, 14, 605, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1013, 3, 14, 15, 608, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1016, 3, 15, 16, 611, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1019, 3, 16, 17, 614, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1022, 0, 593, 1001, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1025, 0, 596, 1004, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1028, 0, 599, 1007, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1031, 0, 602, 1010, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1034, 0, 605, 1013, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1037, 0, 608, 1016, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1040, 0, 611, 1019, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1043, 3, 617, 53, 56, 644, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1046, 3, 620, 56, 59, 647, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1049, 3, 623, 59, 62, 650, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1052, 3, 626, 62, 65, 653, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1055, 3, 629, 65, 68, 656, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1058, 3, 632, 68, 71, 659, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1061, 3, 635, 71, 74, 662, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1064, 3, 638, 74, 77, 665, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1067, 3, 641, 77, 80, 668, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1070, 0, 3, 647, 1049, 86, 92, 674, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1079, 0, 3, 650, 1052, 92, 98, 677, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1088, 0, 3, 653, 1055, 98, 104, 680, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1097, 0, 3, 656, 1058, 104, 110, 683, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1106, 0, 3, 659, 1061, 110, 116, 686, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1115, 0, 3, 662, 1064, 116, 122, 689, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1124, 0, 3, 665, 1067, 122, 128, 692, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1133, 0, 3, 1043, 1046, 671, 1070, 140, 149, 695, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1148, 0, 3, 1046, 1049, 674, 1079, 149, 158, 698, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1163, 0, 3, 1049, 1052, 677, 1088, 158, 167, 701, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_39(buffer, 1178, 0, 3, 1052, 1055, 680, 1097, 167, 176, 704, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1193, 0, 3, 1055, 1058, 683, 1106, 176, 185, 713, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1214, 0, 3, 1058, 1061, 686, 1115, 185, 194, 722, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1229, 0, 3, 1061, 1064, 689, 1124, 194, 203, 725, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1244, 0, 3, 1070, 1079, 698, 1163, 221, 233, 731, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_41(buffer, 1265, 0, 3, 1079, 1088, 701, 1178, 233, 245, 734, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_42(buffer, 1289, 0, 3, 1088, 1097, 704, 1193, 245, 257, 749, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_43(buffer, 1331, 0, 3, 1097, 1106, 713, 1214, 257, 269, 764, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_22(buffer, 1361, 0, 3, 1106, 1115, 722, 1229, 269, 281, 773, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_28(buffer, 1385, 0, 3, 1133, 1148, 728, 1244, 305, 320, 776, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_29(buffer, 1412, 0, 3, 1148, 1163, 731, 1265, 320, 335, 779, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_30(buffer, 1439, 0, 3, 1163, 1178, 734, 1289, 335, 350, 797, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_31(buffer, 1514, 0, 3, 1178, 1193, 749, 1331, 350, 365, 815, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_32(buffer, 1571, 0, 3, 1193, 1214, 764, 1361, 365, 380, 833, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_12(buffer, 1613, 0, 3, 1244, 1265, 779, 1439, 410, 428, 863, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_13(buffer, 1664, 0, 3, 1265, 1289, 797, 1514, 428, 446, 884, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_14(buffer, 1754, 0, 3, 1289, 1331, 815, 1571, 446, 464, 905, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_2(buffer, 1829, 0, 3, 1385, 1412, 842, 1613, 500, 518, 926, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_3(buffer, 1883, 0, 3, 1412, 1439, 863, 1664, 518, 536, 950, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_4(buffer, 1937, 0, 3, 1439, 1514, 884, 1754, 536, 554, 974, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2051, 3, 590, 593, 1001, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2054, 3, 593, 596, 1004, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2057, 3, 596, 599, 1007, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2060, 3, 599, 602, 1010, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2063, 3, 602, 605, 1013, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2066, 3, 605, 608, 1016, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2069, 3, 608, 611, 1019, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2072, 0, 998, 2051, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2075, 0, 1001, 2054, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2078, 0, 1004, 2057, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2081, 0, 1007, 2060, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2084, 0, 1010, 2063, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2087, 0, 1013, 2066, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2090, 0, 1016, 2069, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2093, 3, 1022, 644, 647, 1049, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2096, 3, 1025, 647, 650, 1052, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2099, 3, 1028, 650, 653, 1055, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2102, 3, 1031, 653, 656, 1058, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2105, 3, 1034, 656, 659, 1061, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2108, 3, 1037, 659, 662, 1064, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2111, 3, 1040, 662, 665, 1067, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 2114, 0, 3, 1049, 2096, 671, 674, 1079, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 2126, 0, 3, 1052, 2099, 674, 677, 1088, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 2138, 0, 3, 1055, 2102, 677, 680, 1097, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 2150, 0, 3, 1058, 2105, 680, 683, 1106, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_27(buffer, 2162, 0, 3, 1061, 2108, 683, 686, 1115, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2174, 0, 3, 1064, 2111, 686, 689, 1124, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_53(buffer, 2183, 0, 3, 2093, 2096, 1079, 2126, 695, 698, 1163, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_41(buffer, 2201, 0, 3, 2096, 2099, 1088, 2138, 698, 701, 1178, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_54(buffer, 2222, 0, 3, 2099, 2102, 1097, 2150, 701, 704, 1193, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_55(buffer, 2243, 0, 3, 2102, 2105, 1106, 2162, 704, 713, 1214, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_56(buffer, 2264, 0, 3, 2105, 2108, 1115, 2174, 713, 722, 1229, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_43(buffer, 2282, 0, 3, 2114, 2126, 1163, 2201, 728, 731, 1265, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_44(buffer, 2315, 0, 3, 2126, 2138, 1178, 2222, 731, 734, 1289, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_45(buffer, 2348, 0, 3, 2138, 2150, 1193, 2243, 734, 749, 1331, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_46(buffer, 2387, 0, 3, 2150, 2162, 1214, 2264, 749, 764, 1361, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_26(buffer, 2417, 0, 3, 2183, 2201, 1265, 2315, 776, 779, 1439, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_27(buffer, 2462, 0, 3, 2201, 2222, 1289, 2348, 779, 797, 1514, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_28(buffer, 2554, 0, 3, 2222, 2243, 1331, 2387, 797, 815, 1571, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_11(buffer, 2614, 0, 3, 2282, 2315, 1439, 2462, 842, 863, 1664, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_12(buffer, 2852, 0, 3, 2315, 2348, 1514, 2554, 863, 884, 1754, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_1(buffer, 2989, 0, 3, 2417, 2462, 1664, 2852, 926, 950, 1937, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3292, 3, 2072, 1043, 1046, 2093, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3295, 3, 2075, 1046, 1049, 2096, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3298, 3, 2078, 1049, 1052, 2099, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3301, 3, 2081, 1052, 1055, 2102, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3304, 3, 2084, 1055, 1058, 2105, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3307, 3, 2087, 1058, 1061, 2108, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3310, 3, 2090, 1061, 1064, 2111, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 3313, 0, 3, 2096, 3298, 1070, 1079, 2126, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 3322, 0, 3, 2099, 3301, 1079, 1088, 2138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 3331, 0, 3, 2102, 3304, 1088, 1097, 2150, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_27(buffer, 3340, 0, 3, 2105, 3307, 1097, 1106, 2162, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 3349, 0, 3, 2108, 3310, 1106, 1115, 2174, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_23(buffer, 3358, 0, 3, 3292, 3295, 2114, 3313, 1133, 1148, 2183, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_39(buffer, 3376, 0, 3, 3295, 3298, 2126, 3322, 1148, 1163, 2201, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_39(buffer, 3394, 0, 3, 3298, 3301, 2138, 3331, 1163, 1178, 2222, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_55(buffer, 3412, 0, 3, 3301, 3304, 2150, 3340, 1178, 1193, 2243, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_56(buffer, 3430, 0, 3, 3304, 3307, 2162, 3349, 1193, 1214, 2264, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_41(buffer, 3448, 0, 3, 3313, 3322, 2201, 3394, 1244, 1265, 2315, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_42(buffer, 3481, 0, 3, 3322, 3331, 2222, 3412, 1265, 1289, 2348, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_43(buffer, 3514, 0, 3, 3331, 3340, 2243, 3430, 1289, 1331, 2387, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_23(buffer, 3547, 0, 3, 3358, 3376, 2282, 3448, 1385, 1412, 2417, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_24(buffer, 3601, 0, 3, 3376, 3394, 2315, 3481, 1412, 1439, 2462, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_25(buffer, 3655, 0, 3, 3394, 3412, 2348, 3514, 1439, 1514, 2554, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_10(buffer, 3724, 0, 3, 3448, 3481, 2462, 3655, 1613, 1664, 2852, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 3913, 0, 3, 3547, 3601, 2614, 3724, 1829, 1883, 2989, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4588, 3913, 675, ncols);
        }
    }

    simdtrf::transform_lg(values, nvalues, buffer, 4588, nmax);
}

}  // namespace simdt2ceri
