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


#include "SimdElectronRepulsionRecLI.hpp"

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
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
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
#include "SimdTransformLI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_li_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_li_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(12490, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 14, ncols, fj, mu);

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

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 64, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 67, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 73, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 76, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 79, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 82, 0, 14, 15, 43, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 85, 0, 15, 16, 46, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 88, 0, 16, 17, 49, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 91, 0, 17, 18, 52, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 94, 0, 18, 19, 55, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 97, 0, 19, 20, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 22, 25, 67, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 106, 0, 25, 28, 70, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 112, 0, 28, 31, 73, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 118, 0, 31, 34, 76, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 124, 0, 34, 37, 79, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 130, 0, 37, 40, 82, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 136, 0, 40, 43, 85, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 142, 0, 43, 46, 88, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 148, 0, 46, 49, 91, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 154, 0, 49, 52, 94, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 160, 0, 52, 55, 97, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 166, 0, 61, 64, 100, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 175, 0, 64, 67, 106, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 184, 0, 67, 70, 112, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 193, 0, 70, 73, 118, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 202, 0, 73, 76, 124, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 211, 0, 76, 79, 130, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 220, 0, 79, 82, 136, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 229, 0, 82, 85, 142, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 238, 0, 85, 88, 148, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 247, 0, 88, 91, 154, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 256, 0, 91, 94, 160, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 265, 0, 100, 106, 184, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 277, 0, 106, 112, 193, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 289, 0, 112, 118, 202, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 301, 0, 118, 124, 211, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 313, 0, 124, 130, 220, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 325, 0, 130, 136, 229, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 337, 0, 136, 142, 238, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 349, 0, 142, 148, 247, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 361, 0, 148, 154, 256, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 373, 0, 166, 175, 265, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 388, 0, 175, 184, 277, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 403, 0, 184, 193, 289, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 418, 0, 193, 202, 301, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 433, 0, 202, 211, 313, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 448, 0, 211, 220, 325, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 463, 0, 220, 229, 337, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 478, 0, 229, 238, 349, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 493, 0, 238, 247, 361, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 508, 0, 265, 277, 403, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 526, 0, 277, 289, 418, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 544, 0, 289, 301, 433, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 562, 0, 301, 313, 448, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 580, 0, 313, 325, 463, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 598, 0, 325, 337, 478, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 616, 0, 337, 349, 493, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 634, 0, 373, 388, 508, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 652, 0, 388, 403, 526, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 670, 0, 403, 418, 544, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 688, 0, 418, 433, 562, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 706, 0, 433, 448, 580, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 724, 0, 448, 463, 598, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 742, 0, 463, 478, 616, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 760, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 763, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 766, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 769, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 772, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 775, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 778, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 781, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 784, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 787, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 790, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 793, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 796, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 799, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 802, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 805, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 808, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 811, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 814, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 817, 3, 17, 49, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 820, 3, 18, 52, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 823, 3, 19, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 826, 3, 25, 67, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 829, 3, 28, 70, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 832, 3, 31, 73, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 835, 3, 34, 76, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 838, 3, 37, 79, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 841, 3, 40, 82, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 844, 3, 43, 85, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 847, 3, 46, 88, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 850, 3, 49, 91, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 853, 3, 52, 94, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 856, 3, 55, 97, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 859, 3, 67, 106, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 862, 3, 70, 112, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 865, 3, 73, 118, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 868, 3, 76, 124, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 871, 3, 79, 130, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 874, 3, 82, 136, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 877, 3, 85, 142, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 880, 3, 88, 148, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 883, 3, 91, 154, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 886, 3, 94, 160, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 889, 3, 106, 184, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 892, 3, 112, 193, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 895, 3, 118, 202, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 898, 3, 124, 211, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 907, 3, 130, 220, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 916, 3, 136, 229, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 925, 3, 142, 238, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 934, 3, 148, 247, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 937, 3, 154, 256, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 940, 3, 184, 277, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 943, 3, 193, 289, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 946, 3, 202, 301, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 961, 3, 211, 313, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 976, 3, 220, 325, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 991, 3, 229, 337, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 1006, 3, 238, 349, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1015, 3, 247, 361, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 1018, 3, 277, 403, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1021, 3, 289, 418, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1039, 3, 301, 433, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1057, 3, 313, 448, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1075, 3, 325, 463, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1093, 3, 337, 478, ncols, p);

            compute_prim_ip_electron_repulsion_17(buffer, 1111, 3, 349, 493, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1120, 3, 403, 526, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1141, 3, 418, 544, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1162, 3, 433, 562, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1183, 3, 448, 580, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1204, 3, 463, 598, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1225, 3, 478, 616, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1246, 3, 526, 670, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1270, 3, 544, 688, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1294, 3, 562, 706, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1318, 3, 580, 724, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1342, 3, 598, 742, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1366, 3, 9, 10, 763, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1369, 3, 10, 11, 766, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1372, 3, 11, 12, 769, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1375, 3, 12, 13, 772, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1378, 3, 13, 14, 775, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1381, 3, 14, 15, 778, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1384, 3, 15, 16, 781, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1387, 3, 16, 17, 784, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1390, 3, 17, 18, 787, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1393, 3, 18, 19, 790, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1396, 0, 763, 1369, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1399, 0, 766, 1372, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1402, 0, 769, 1375, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1405, 0, 772, 1378, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1408, 0, 775, 1381, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1411, 0, 778, 1384, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1414, 0, 781, 1387, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1417, 0, 784, 1390, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1420, 0, 787, 1393, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1423, 3, 793, 61, 64, 826, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1426, 3, 796, 64, 67, 829, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1429, 3, 799, 67, 70, 832, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1432, 3, 802, 70, 73, 835, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1435, 3, 805, 73, 76, 838, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1438, 3, 808, 76, 79, 841, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1441, 3, 811, 79, 82, 844, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1444, 3, 814, 82, 85, 847, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1447, 3, 817, 85, 88, 850, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1450, 3, 820, 88, 91, 853, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1453, 3, 823, 91, 94, 856, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1456, 0, 3, 829, 1429, 100, 106, 862, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1465, 0, 3, 832, 1432, 106, 112, 865, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1474, 0, 3, 835, 1435, 112, 118, 868, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1483, 0, 3, 838, 1438, 118, 124, 871, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1492, 0, 3, 841, 1441, 124, 130, 874, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1501, 0, 3, 844, 1444, 130, 136, 877, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1510, 0, 3, 847, 1447, 136, 142, 880, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1519, 0, 3, 850, 1450, 142, 148, 883, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1528, 0, 3, 853, 1453, 148, 154, 886, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1537, 0, 3, 1423, 1426, 859, 1456, 166, 175, 889, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1552, 0, 3, 1426, 1429, 862, 1465, 175, 184, 892, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1567, 0, 3, 1429, 1432, 865, 1474, 184, 193, 895, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_39(buffer, 1582, 0, 3, 1432, 1435, 868, 1483, 193, 202, 898, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1597, 0, 3, 1435, 1438, 871, 1492, 202, 211, 907, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1618, 0, 3, 1438, 1441, 874, 1501, 211, 220, 916, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1639, 0, 3, 1441, 1444, 877, 1510, 220, 229, 925, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1660, 0, 3, 1444, 1447, 880, 1519, 229, 238, 934, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1675, 0, 3, 1447, 1450, 883, 1528, 238, 247, 937, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1690, 0, 3, 1456, 1465, 892, 1567, 265, 277, 943, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1711, 0, 3, 1465, 1474, 895, 1582, 277, 289, 946, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1732, 0, 3, 1474, 1483, 898, 1597, 289, 301, 961, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1765, 0, 3, 1483, 1492, 907, 1618, 301, 313, 976, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1798, 0, 3, 1492, 1501, 916, 1639, 313, 325, 991, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_45(buffer, 1831, 0, 3, 1501, 1510, 925, 1660, 325, 337, 1006, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1858, 0, 3, 1510, 1519, 934, 1675, 337, 349, 1015, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_28(buffer, 1879, 0, 3, 1537, 1552, 940, 1690, 373, 388, 1018, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_33(buffer, 1906, 0, 3, 1552, 1567, 943, 1711, 388, 403, 1021, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 1933, 0, 3, 1567, 1582, 946, 1732, 403, 418, 1039, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_37(buffer, 1975, 0, 3, 1582, 1597, 961, 1765, 418, 433, 1057, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2017, 0, 3, 1597, 1618, 976, 1798, 433, 448, 1075, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_35(buffer, 2059, 0, 3, 1618, 1639, 991, 1831, 448, 463, 1093, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_36(buffer, 2101, 0, 3, 1639, 1660, 1006, 1858, 463, 478, 1111, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_19(buffer, 2134, 0, 3, 1690, 1711, 1021, 1933, 508, 526, 1141, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_20(buffer, 2185, 0, 3, 1711, 1732, 1039, 1975, 526, 544, 1162, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2236, 0, 3, 1732, 1765, 1057, 2017, 544, 562, 1183, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2287, 0, 3, 1765, 1798, 1075, 2059, 562, 580, 1204, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_18(buffer, 2338, 0, 3, 1798, 1831, 1093, 2101, 580, 598, 1225, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_2(buffer, 2389, 0, 3, 1879, 1906, 1120, 2134, 634, 652, 1246, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_5(buffer, 2443, 0, 3, 1906, 1933, 1141, 2185, 652, 670, 1270, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2497, 0, 3, 1933, 1975, 1162, 2236, 670, 688, 1294, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2551, 0, 3, 1975, 2017, 1183, 2287, 688, 706, 1318, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2605, 0, 3, 2017, 2059, 1204, 2338, 706, 724, 1342, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2659, 3, 760, 763, 1369, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2662, 3, 763, 766, 1372, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2665, 3, 766, 769, 1375, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2668, 3, 769, 772, 1378, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2671, 3, 772, 775, 1381, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2674, 3, 775, 778, 1384, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2677, 3, 778, 781, 1387, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2680, 3, 781, 784, 1390, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2683, 3, 784, 787, 1393, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2686, 0, 1366, 2659, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2689, 0, 1369, 2662, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2692, 0, 1372, 2665, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2695, 0, 1375, 2668, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2698, 0, 1378, 2671, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2701, 0, 1381, 2674, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2704, 0, 1384, 2677, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2707, 0, 1387, 2680, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2710, 0, 1390, 2683, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2713, 3, 1396, 826, 829, 1429, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2716, 3, 1399, 829, 832, 1432, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2722, 3, 1402, 832, 835, 1435, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2728, 3, 1405, 835, 838, 1438, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2734, 3, 1408, 838, 841, 1441, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2740, 3, 1411, 841, 844, 1444, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2746, 3, 1414, 844, 847, 1447, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2752, 3, 1417, 847, 850, 1450, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2758, 3, 1420, 850, 853, 1453, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 2761, 0, 3, 1429, 2716, 859, 862, 1465, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 2773, 0, 3, 1432, 2722, 862, 865, 1474, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2785, 0, 3, 1435, 2728, 865, 868, 1483, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2800, 0, 3, 1438, 2734, 868, 871, 1492, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2815, 0, 3, 1441, 2740, 871, 874, 1501, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2830, 0, 3, 1444, 2746, 874, 877, 1510, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2845, 0, 3, 1447, 2752, 877, 880, 1519, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2860, 0, 3, 1450, 2758, 880, 883, 1528, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_57(buffer, 2869, 0, 3, 2713, 2716, 1465, 2773, 889, 892, 1567, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_58(buffer, 2887, 0, 3, 2716, 2722, 1474, 2785, 892, 895, 1582, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_59(buffer, 2911, 0, 3, 2722, 2728, 1483, 2800, 895, 898, 1597, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 2935, 0, 3, 2728, 2734, 1492, 2815, 898, 907, 1618, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 2965, 0, 3, 2734, 2740, 1501, 2830, 907, 916, 1639, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_61(buffer, 2995, 0, 3, 2740, 2746, 1510, 2845, 916, 925, 1660, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_62(buffer, 3019, 0, 3, 2746, 2752, 1519, 2860, 925, 934, 1675, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_52(buffer, 3037, 0, 3, 2761, 2773, 1567, 2887, 940, 943, 1711, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_53(buffer, 3067, 0, 3, 2773, 2785, 1582, 2911, 943, 946, 1732, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_54(buffer, 3097, 0, 3, 2785, 2800, 1597, 2935, 946, 961, 1765, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 3151, 0, 3, 2800, 2815, 1618, 2965, 961, 976, 1798, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_56(buffer, 3202, 0, 3, 2815, 2830, 1639, 2995, 976, 991, 1831, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_57(buffer, 3247, 0, 3, 2830, 2845, 1660, 3019, 991, 1006, 1858, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_34(buffer, 3274, 0, 3, 2869, 2887, 1711, 3067, 1018, 1021, 1933, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_35(buffer, 3310, 0, 3, 2887, 2911, 1732, 3097, 1021, 1039, 1975, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_36(buffer, 3391, 0, 3, 2911, 2935, 1765, 3151, 1039, 1057, 2017, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_37(buffer, 3478, 0, 3, 2935, 2965, 1798, 3202, 1057, 1075, 2059, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_38(buffer, 3559, 0, 3, 2965, 2995, 1831, 3247, 1075, 1093, 2101, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_16(buffer, 3610, 0, 3, 3037, 3067, 1933, 3310, 1120, 1141, 2185, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_17(buffer, 3706, 0, 3, 3067, 3097, 1975, 3391, 1141, 1162, 2236, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_18(buffer, 3805, 0, 3, 3097, 3151, 2017, 3478, 1162, 1183, 2287, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_19(buffer, 3907, 0, 3, 3151, 3202, 2059, 3559, 1183, 1204, 2338, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_5(buffer, 4003, 0, 3, 3274, 3310, 2185, 3706, 1246, 1270, 2497, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_6(buffer, 4111, 0, 3, 3310, 3391, 2236, 3805, 1270, 1294, 2551, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_7(buffer, 4219, 0, 3, 3391, 3478, 2287, 3907, 1294, 1318, 2605, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4327, 3, 1366, 1369, 2662, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4330, 3, 1369, 1372, 2665, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4333, 3, 1372, 1375, 2668, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4336, 3, 1375, 1378, 2671, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4339, 3, 1378, 1381, 2674, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4342, 3, 1381, 1384, 2677, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4345, 3, 1384, 1387, 2680, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4348, 3, 1387, 1390, 2683, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4351, 0, 2662, 4330, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4354, 0, 2665, 4333, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4357, 0, 2668, 4336, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4360, 0, 2671, 4339, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4363, 0, 2674, 4342, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4366, 0, 2677, 4345, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4369, 0, 2680, 4348, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 4372, 3, 2686, 1423, 1426, 2713, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 4375, 3, 2689, 1426, 1429, 2716, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4378, 3, 2692, 1429, 1432, 2722, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4387, 3, 2695, 1432, 1435, 2728, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4396, 3, 2698, 1435, 1438, 2734, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4405, 3, 2701, 1438, 1441, 2740, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4414, 3, 2704, 1441, 1444, 2746, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4423, 3, 2707, 1444, 1447, 2752, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 4432, 3, 2710, 1447, 1450, 2758, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_34(buffer, 4435, 0, 3, 2716, 4378, 1456, 1465, 2773, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_55(buffer, 4444, 0, 3, 2722, 4387, 1465, 1474, 2785, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 4459, 0, 3, 2728, 4396, 1474, 1483, 2800, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 4480, 0, 3, 2734, 4405, 1483, 1492, 2815, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 4501, 0, 3, 2740, 4414, 1492, 1501, 2830, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 4522, 0, 3, 2746, 4423, 1501, 1510, 2845, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_54(buffer, 4543, 0, 3, 2752, 4432, 1510, 1519, 2860, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_29(buffer, 4552, 0, 3, 4372, 4375, 2761, 4435, 1537, 1552, 2869, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_62(buffer, 4573, 0, 3, 4375, 4378, 2773, 4444, 1552, 1567, 2887, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_63(buffer, 4594, 0, 3, 4378, 4387, 2785, 4459, 1567, 1582, 2911, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_64(buffer, 4627, 0, 3, 4387, 4396, 2800, 4480, 1582, 1597, 2935, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 4660, 0, 3, 4396, 4405, 2815, 4501, 1597, 1618, 2965, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_66(buffer, 4699, 0, 3, 4405, 4414, 2830, 4522, 1618, 1639, 2995, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_67(buffer, 4732, 0, 3, 4414, 4423, 2845, 4543, 1639, 1660, 3019, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_48(buffer, 4753, 0, 3, 4435, 4444, 2887, 4594, 1690, 1711, 3067, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_49(buffer, 4787, 0, 3, 4444, 4459, 2911, 4627, 1711, 1732, 3097, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_50(buffer, 4830, 0, 3, 4459, 4480, 2935, 4660, 1732, 1765, 3151, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_51(buffer, 4909, 0, 3, 4480, 4501, 2965, 4699, 1765, 1798, 3202, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_52(buffer, 4976, 0, 3, 4501, 4522, 2995, 4732, 1798, 1831, 3247, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_29(buffer, 5013, 0, 3, 4552, 4573, 3037, 4753, 1879, 1906, 3274, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_30(buffer, 5061, 0, 3, 4573, 4594, 3067, 4787, 1906, 1933, 3310, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_31(buffer, 5109, 0, 3, 4594, 4627, 3097, 4830, 1933, 1975, 3391, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_32(buffer, 5274, 0, 3, 4627, 4660, 3151, 4909, 1975, 2017, 3478, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_33(buffer, 5415, 0, 3, 4660, 4699, 3202, 4976, 2017, 2059, 3559, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_13(buffer, 5496, 0, 3, 4753, 4787, 3310, 5109, 2134, 2185, 3706, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_14(buffer, 5643, 0, 3, 4787, 4830, 3391, 5274, 2185, 2236, 3805, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_15(buffer, 5846, 0, 3, 4830, 4909, 3478, 5415, 2236, 2287, 3907, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_2(buffer, 6017, 0, 3, 5013, 5061, 3610, 5496, 2389, 2443, 4003, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_3(buffer, 6185, 0, 3, 5061, 5109, 3706, 5643, 2443, 2497, 4111, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_4(buffer, 6353, 0, 3, 5109, 5274, 3805, 5846, 2497, 2551, 4219, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6581, 3, 2659, 2662, 4330, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6584, 3, 2662, 2665, 4333, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6587, 3, 2665, 2668, 4336, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6590, 3, 2668, 2671, 4339, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6593, 3, 2671, 2674, 4342, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6596, 3, 2674, 2677, 4345, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6599, 3, 2677, 2680, 4348, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6602, 0, 4327, 6581, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6605, 0, 4330, 6584, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6608, 0, 4333, 6587, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6611, 0, 4336, 6590, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6614, 0, 4339, 6593, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6617, 0, 4342, 6596, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6620, 0, 4345, 6599, ncols, p);

            compute_prim_dh_electron_repulsion_32(buffer, 6623, 3, 4351, 2713, 2716, 4378, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 6626, 3, 4354, 2716, 2722, 4387, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 6635, 3, 4357, 2722, 2728, 4396, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 6644, 3, 4360, 2728, 2734, 4405, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 6653, 3, 4363, 2734, 2740, 4414, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 6662, 3, 4366, 2740, 2746, 4423, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 6671, 3, 4369, 2746, 2752, 4432, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_30(buffer, 6674, 0, 3, 4378, 6626, 2761, 2773, 4444, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_57(buffer, 6692, 0, 3, 4387, 6635, 2773, 2785, 4459, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_47(buffer, 6710, 0, 3, 4396, 6644, 2785, 2800, 4480, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_47(buffer, 6734, 0, 3, 4405, 6653, 2800, 2815, 4501, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_47(buffer, 6758, 0, 3, 4414, 6662, 2815, 2830, 4522, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_58(buffer, 6782, 0, 3, 4423, 6671, 2830, 2845, 4543, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_52(buffer, 6791, 0, 3, 6623, 6626, 4444, 6692, 2869, 2887, 4594, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_53(buffer, 6815, 0, 3, 6626, 6635, 4459, 6710, 2887, 2911, 4627, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_54(buffer, 6854, 0, 3, 6635, 6644, 4480, 6734, 2911, 2935, 4660, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_55(buffer, 6893, 0, 3, 6644, 6653, 4501, 6758, 2935, 2965, 4699, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_56(buffer, 6932, 0, 3, 6653, 6662, 4522, 6782, 2965, 2995, 4732, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_36(buffer, 6956, 0, 3, 6674, 6692, 4594, 6815, 3037, 3067, 4787, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_37(buffer, 7013, 0, 3, 6692, 6710, 4627, 6854, 3067, 3097, 4830, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_38(buffer, 7070, 0, 3, 6710, 6734, 4660, 6893, 3097, 3151, 4909, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_39(buffer, 7151, 0, 3, 6734, 6758, 4699, 6932, 3151, 3202, 4976, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_21(buffer, 7199, 0, 3, 6791, 6815, 4787, 7013, 3274, 3310, 5109, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_22(buffer, 7280, 0, 3, 6815, 6854, 4830, 7070, 3310, 3391, 5274, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_23(buffer, 7483, 0, 3, 6854, 6893, 4909, 7151, 3391, 3478, 5415, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_9(buffer, 7600, 0, 3, 6956, 7013, 5109, 7280, 3610, 3706, 5643, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_10(buffer, 8081, 0, 3, 7013, 7070, 5274, 7483, 3706, 3805, 5846, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_1(buffer, 8377, 0, 3, 7199, 7280, 5643, 8081, 4003, 4111, 6353, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_9(buffer, 8974, 3, 6602, 4372, 4375, 6623, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_25(buffer, 8977, 3, 6605, 4375, 4378, 6626, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 8980, 3, 6608, 4378, 4387, 6635, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 8983, 3, 6611, 4387, 4396, 6644, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 8986, 3, 6614, 4396, 4405, 6653, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 8989, 3, 6617, 4405, 4414, 6662, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 8992, 3, 6620, 4414, 4423, 6671, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_20(buffer, 8995, 0, 3, 6626, 8980, 4435, 4444, 6692, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_44(buffer, 9004, 0, 3, 6635, 8983, 4444, 4459, 6710, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_32(buffer, 9013, 0, 3, 6644, 8986, 4459, 4480, 6734, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_32(buffer, 9022, 0, 3, 6653, 8989, 4480, 4501, 6758, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_45(buffer, 9031, 0, 3, 6662, 8992, 4501, 4522, 6782, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_11(buffer, 9040, 0, 3, 8974, 8977, 6674, 8995, 4552, 4573, 6791, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_34(buffer, 9064, 0, 3, 8977, 8980, 6692, 9004, 4573, 4594, 6815, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_35(buffer, 9088, 0, 3, 8980, 8983, 6710, 9013, 4594, 4627, 6854, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_36(buffer, 9112, 0, 3, 8983, 8986, 6734, 9022, 4627, 4660, 6893, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_37(buffer, 9136, 0, 3, 8986, 8989, 6758, 9031, 4660, 4699, 6932, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_22(buffer, 9160, 0, 3, 8995, 9004, 6815, 9088, 4753, 4787, 7013, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_23(buffer, 9211, 0, 3, 9004, 9013, 6854, 9112, 4787, 4830, 7070, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_24(buffer, 9262, 0, 3, 9013, 9022, 6893, 9136, 4830, 4909, 7151, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_11(buffer, 9313, 0, 3, 9040, 9064, 6956, 9160, 5013, 5061, 7199, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_12(buffer, 9403, 0, 3, 9064, 9088, 7013, 9211, 5061, 5109, 7280, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_13(buffer, 9493, 0, 3, 9088, 9112, 7070, 9262, 5109, 5274, 7483, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_5(buffer, 9619, 0, 3, 9160, 9211, 7280, 9493, 5496, 5643, 8081, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 9970, 0, 3, 9313, 9403, 7600, 9619, 6017, 6185, 8377, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 11230, 9970, 1260, ncols);
        }
    }

    simdtrf::transform_li(values, nvalues, buffer, 11230, nmax);
}

}  // namespace simdt2ceri
