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


#include "SimdElectronRepulsionRecLK.hpp"

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
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformLK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(18119, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 64, 0, 7, 8, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 67, 0, 8, 9, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 70, 0, 9, 10, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 73, 0, 10, 11, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 76, 0, 11, 12, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 79, 0, 12, 13, 40, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 82, 0, 13, 14, 43, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 85, 0, 14, 15, 46, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 88, 0, 15, 16, 49, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 91, 0, 16, 17, 52, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 94, 0, 17, 18, 55, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 97, 0, 18, 19, 58, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 100, 0, 19, 20, 61, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 103, 0, 22, 25, 67, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 109, 0, 25, 28, 70, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 115, 0, 28, 31, 73, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 121, 0, 31, 34, 76, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 127, 0, 34, 37, 79, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 133, 0, 37, 40, 82, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 139, 0, 40, 43, 85, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 145, 0, 43, 46, 88, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 151, 0, 46, 49, 91, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 157, 0, 49, 52, 94, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 163, 0, 52, 55, 97, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 169, 0, 55, 58, 100, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 175, 0, 64, 67, 109, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 184, 0, 67, 70, 115, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 193, 0, 70, 73, 121, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 202, 0, 73, 76, 127, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 211, 0, 76, 79, 133, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 220, 0, 79, 82, 139, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 229, 0, 82, 85, 145, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 238, 0, 85, 88, 151, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 247, 0, 88, 91, 157, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 256, 0, 91, 94, 163, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 265, 0, 94, 97, 169, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 274, 0, 103, 109, 184, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 286, 0, 109, 115, 193, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 298, 0, 115, 121, 202, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 310, 0, 121, 127, 211, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 322, 0, 127, 133, 220, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 334, 0, 133, 139, 229, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 346, 0, 139, 145, 238, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 358, 0, 145, 151, 247, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 370, 0, 151, 157, 256, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 382, 0, 157, 163, 265, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 394, 0, 175, 184, 286, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 409, 0, 184, 193, 298, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 424, 0, 193, 202, 310, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 439, 0, 202, 211, 322, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 454, 0, 211, 220, 334, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 469, 0, 220, 229, 346, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 484, 0, 229, 238, 358, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 499, 0, 238, 247, 370, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 514, 0, 247, 256, 382, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 529, 0, 274, 286, 409, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 544, 0, 286, 298, 424, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 562, 0, 298, 310, 439, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 580, 0, 310, 322, 454, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 598, 0, 322, 334, 469, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 616, 0, 334, 346, 484, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 634, 0, 346, 358, 499, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 652, 0, 358, 370, 514, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 670, 0, 394, 409, 544, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 688, 0, 409, 424, 562, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 706, 0, 424, 439, 580, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 724, 0, 439, 454, 598, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 742, 0, 454, 469, 616, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 760, 0, 469, 484, 634, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 778, 0, 484, 499, 652, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 796, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 799, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 802, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 805, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 808, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 811, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 814, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 817, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 820, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 823, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 826, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 829, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 832, 3, 9, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 835, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 838, 3, 11, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 841, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 844, 3, 13, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 847, 3, 14, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 850, 3, 15, 46, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 853, 3, 16, 49, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 856, 3, 17, 52, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 859, 3, 18, 55, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 862, 3, 19, 58, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 865, 3, 22, 64, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 868, 3, 25, 67, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 871, 3, 28, 70, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 874, 3, 31, 73, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 877, 3, 34, 76, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 880, 3, 37, 79, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 883, 3, 40, 82, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 886, 3, 43, 85, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 889, 3, 46, 88, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 892, 3, 49, 91, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 895, 3, 52, 94, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 898, 3, 55, 97, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 901, 3, 58, 100, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 904, 3, 67, 109, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 907, 3, 70, 115, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 910, 3, 73, 121, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 913, 3, 76, 127, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 916, 3, 79, 133, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 919, 3, 82, 139, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 922, 3, 85, 145, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 925, 3, 88, 151, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 928, 3, 91, 157, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 931, 3, 94, 163, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 934, 3, 97, 169, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 937, 3, 103, 175, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 940, 3, 109, 184, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 943, 3, 115, 193, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 946, 3, 121, 202, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 949, 3, 127, 211, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 958, 3, 133, 220, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 967, 3, 139, 229, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 976, 3, 145, 238, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 985, 3, 151, 247, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 994, 3, 157, 256, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 997, 3, 163, 265, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1000, 3, 184, 286, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1003, 3, 193, 298, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1006, 3, 202, 310, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1021, 3, 211, 322, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1036, 3, 220, 334, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1051, 3, 229, 346, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1066, 3, 238, 358, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 1081, 3, 247, 370, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1090, 3, 256, 382, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 1093, 3, 274, 394, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 1096, 3, 286, 409, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1099, 3, 298, 424, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1117, 3, 310, 439, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1135, 3, 322, 454, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1153, 3, 334, 469, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1171, 3, 346, 484, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1189, 3, 358, 499, ncols, p);

            compute_prim_ip_electron_repulsion_17(buffer, 1207, 3, 370, 514, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1216, 3, 409, 544, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1237, 3, 424, 562, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1258, 3, 439, 580, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1279, 3, 454, 598, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1300, 3, 469, 616, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1321, 3, 484, 634, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1342, 3, 499, 652, ncols, p);

            compute_prim_lp_electron_repulsion_2(buffer, 1363, 3, 529, 670, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1387, 3, 544, 688, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1411, 3, 562, 706, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1435, 3, 580, 724, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1459, 3, 598, 742, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1483, 3, 616, 760, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1507, 3, 634, 778, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1531, 3, 8, 9, 799, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1534, 3, 9, 10, 802, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1537, 3, 10, 11, 805, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1540, 3, 11, 12, 808, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1543, 3, 12, 13, 811, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1546, 3, 13, 14, 814, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1549, 3, 14, 15, 817, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1552, 3, 15, 16, 820, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1555, 3, 16, 17, 823, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1558, 3, 17, 18, 826, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1561, 3, 18, 19, 829, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1564, 0, 796, 1531, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1567, 0, 799, 1534, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1570, 0, 802, 1537, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1573, 0, 805, 1540, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1576, 0, 808, 1543, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1579, 0, 811, 1546, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1582, 0, 814, 1549, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1585, 0, 817, 1552, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1588, 0, 820, 1555, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1591, 0, 823, 1558, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1594, 0, 826, 1561, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1597, 3, 832, 64, 67, 871, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1600, 3, 835, 67, 70, 874, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1603, 3, 838, 70, 73, 877, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1606, 3, 841, 73, 76, 880, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1609, 3, 844, 76, 79, 883, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1612, 3, 847, 79, 82, 886, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1615, 3, 850, 82, 85, 889, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1618, 3, 853, 85, 88, 892, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1621, 3, 856, 88, 91, 895, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1624, 3, 859, 91, 94, 898, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1627, 3, 862, 94, 97, 901, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1630, 0, 3, 871, 1600, 103, 109, 907, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1639, 0, 3, 874, 1603, 109, 115, 910, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1648, 0, 3, 877, 1606, 115, 121, 913, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1657, 0, 3, 880, 1609, 121, 127, 916, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1666, 0, 3, 883, 1612, 127, 133, 919, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1675, 0, 3, 886, 1615, 133, 139, 922, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1684, 0, 3, 889, 1618, 139, 145, 925, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1693, 0, 3, 892, 1621, 145, 151, 928, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1702, 0, 3, 895, 1624, 151, 157, 931, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1711, 0, 3, 898, 1627, 157, 163, 934, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1720, 0, 3, 1597, 1600, 907, 1639, 175, 184, 943, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1735, 0, 3, 1600, 1603, 910, 1648, 184, 193, 946, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_39(buffer, 1750, 0, 3, 1603, 1606, 913, 1657, 193, 202, 949, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1765, 0, 3, 1606, 1609, 916, 1666, 202, 211, 958, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1786, 0, 3, 1609, 1612, 919, 1675, 211, 220, 967, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1807, 0, 3, 1612, 1615, 922, 1684, 220, 229, 976, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 1828, 0, 3, 1615, 1618, 925, 1693, 229, 238, 985, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1849, 0, 3, 1618, 1621, 928, 1702, 238, 247, 994, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1864, 0, 3, 1621, 1624, 931, 1711, 247, 256, 997, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1879, 0, 3, 1630, 1639, 943, 1735, 274, 286, 1003, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1900, 0, 3, 1639, 1648, 946, 1750, 286, 298, 1006, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1921, 0, 3, 1648, 1657, 949, 1765, 298, 310, 1021, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1954, 0, 3, 1657, 1666, 958, 1786, 310, 322, 1036, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 1987, 0, 3, 1666, 1675, 967, 1807, 322, 334, 1051, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2020, 0, 3, 1675, 1684, 976, 1828, 334, 346, 1066, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_45(buffer, 2053, 0, 3, 1684, 1693, 985, 1849, 346, 358, 1081, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 2080, 0, 3, 1693, 1702, 994, 1864, 358, 370, 1090, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_33(buffer, 2101, 0, 3, 1720, 1735, 1003, 1900, 394, 409, 1099, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 2128, 0, 3, 1735, 1750, 1006, 1921, 409, 424, 1117, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_37(buffer, 2170, 0, 3, 1750, 1765, 1021, 1954, 424, 439, 1135, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2212, 0, 3, 1765, 1786, 1036, 1987, 439, 454, 1153, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2254, 0, 3, 1786, 1807, 1051, 2020, 454, 469, 1171, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_35(buffer, 2296, 0, 3, 1807, 1828, 1066, 2053, 469, 484, 1189, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_36(buffer, 2338, 0, 3, 1828, 1849, 1081, 2080, 484, 499, 1207, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_15(buffer, 2371, 0, 3, 1879, 1900, 1099, 2128, 529, 544, 1237, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_20(buffer, 2419, 0, 3, 1900, 1921, 1117, 2170, 544, 562, 1258, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2470, 0, 3, 1921, 1954, 1135, 2212, 562, 580, 1279, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2521, 0, 3, 1954, 1987, 1153, 2254, 580, 598, 1300, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2572, 0, 3, 1987, 2020, 1171, 2296, 598, 616, 1321, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_18(buffer, 2623, 0, 3, 2020, 2053, 1189, 2338, 616, 634, 1342, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_5(buffer, 2674, 0, 3, 2101, 2128, 1237, 2419, 670, 688, 1411, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2728, 0, 3, 2128, 2170, 1258, 2470, 688, 706, 1435, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2782, 0, 3, 2170, 2212, 1279, 2521, 706, 724, 1459, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2836, 0, 3, 2212, 2254, 1300, 2572, 724, 742, 1483, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 2890, 0, 3, 2254, 2296, 1321, 2623, 742, 760, 1507, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2944, 3, 796, 799, 1534, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2947, 3, 799, 802, 1537, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2950, 3, 802, 805, 1540, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2953, 3, 805, 808, 1543, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2956, 3, 808, 811, 1546, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2959, 3, 811, 814, 1549, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2962, 3, 814, 817, 1552, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2965, 3, 817, 820, 1555, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2968, 3, 820, 823, 1558, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2971, 3, 823, 826, 1561, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2974, 0, 1534, 2947, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2977, 0, 1537, 2950, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2980, 0, 1540, 2953, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2983, 0, 1543, 2956, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2986, 0, 1546, 2959, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2989, 0, 1549, 2962, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2992, 0, 1552, 2965, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2995, 0, 1555, 2968, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2998, 0, 1558, 2971, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 3001, 3, 1564, 865, 868, 1597, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 3004, 3, 1567, 868, 871, 1600, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3007, 3, 1570, 871, 874, 1603, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3013, 3, 1573, 874, 877, 1606, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3019, 3, 1576, 877, 880, 1609, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3025, 3, 1579, 880, 883, 1612, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3031, 3, 1582, 883, 886, 1615, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3037, 3, 1585, 886, 889, 1618, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3043, 3, 1588, 889, 892, 1621, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3049, 3, 1591, 892, 895, 1624, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 3055, 3, 1594, 895, 898, 1627, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_31(buffer, 3058, 0, 3, 1600, 3007, 904, 907, 1639, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 3067, 0, 3, 1603, 3013, 907, 910, 1648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3079, 0, 3, 1606, 3019, 910, 913, 1657, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3094, 0, 3, 1609, 3025, 913, 916, 1666, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3109, 0, 3, 1612, 3031, 916, 919, 1675, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3124, 0, 3, 1615, 3037, 919, 922, 1684, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3139, 0, 3, 1618, 3043, 922, 925, 1693, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3154, 0, 3, 1621, 3049, 925, 928, 1702, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 3169, 0, 3, 1624, 3055, 928, 931, 1711, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_29(buffer, 3178, 0, 3, 3001, 3004, 1630, 3058, 937, 940, 1720, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_57(buffer, 3196, 0, 3, 3004, 3007, 1639, 3067, 940, 943, 1735, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_58(buffer, 3214, 0, 3, 3007, 3013, 1648, 3079, 943, 946, 1750, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_59(buffer, 3238, 0, 3, 3013, 3019, 1657, 3094, 946, 949, 1765, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3262, 0, 3, 3019, 3025, 1666, 3109, 949, 958, 1786, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3292, 0, 3, 3025, 3031, 1675, 3124, 958, 967, 1807, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3322, 0, 3, 3031, 3037, 1684, 3139, 967, 976, 1828, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_61(buffer, 3352, 0, 3, 3037, 3043, 1693, 3154, 976, 985, 1849, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_62(buffer, 3376, 0, 3, 3043, 3049, 1702, 3169, 985, 994, 1864, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_47(buffer, 3394, 0, 3, 3058, 3067, 1735, 3214, 1000, 1003, 1900, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_53(buffer, 3421, 0, 3, 3067, 3079, 1750, 3238, 1003, 1006, 1921, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 3451, 0, 3, 3079, 3094, 1765, 3262, 1006, 1021, 1954, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 3502, 0, 3, 3094, 3109, 1786, 3292, 1021, 1036, 1987, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 3553, 0, 3, 3109, 3124, 1807, 3322, 1036, 1051, 2020, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_56(buffer, 3604, 0, 3, 3124, 3139, 1828, 3352, 1051, 1066, 2053, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_57(buffer, 3649, 0, 3, 3139, 3154, 1849, 3376, 1066, 1081, 2080, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_29(buffer, 3676, 0, 3, 3178, 3196, 1879, 3394, 1093, 1096, 2101, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_34(buffer, 3712, 0, 3, 3196, 3214, 1900, 3421, 1096, 1099, 2128, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_39(buffer, 3748, 0, 3, 3214, 3238, 1921, 3451, 1099, 1117, 2170, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_40(buffer, 3829, 0, 3, 3238, 3262, 1954, 3502, 1117, 1135, 2212, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_41(buffer, 3910, 0, 3, 3262, 3292, 1987, 3553, 1135, 1153, 2254, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_37(buffer, 3991, 0, 3, 3292, 3322, 2020, 3604, 1153, 1171, 2296, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_38(buffer, 4072, 0, 3, 3322, 3352, 2053, 3649, 1171, 1189, 2338, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_20(buffer, 4123, 0, 3, 3394, 3421, 2128, 3748, 1216, 1237, 2419, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_21(buffer, 4219, 0, 3, 3421, 3451, 2170, 3829, 1237, 1258, 2470, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_22(buffer, 4318, 0, 3, 3451, 3502, 2212, 3910, 1258, 1279, 2521, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_22(buffer, 4417, 0, 3, 3502, 3553, 2254, 3991, 1279, 1300, 2572, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_19(buffer, 4516, 0, 3, 3553, 3604, 2296, 4072, 1300, 1321, 2623, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_2(buffer, 4612, 0, 3, 3676, 3712, 2371, 4123, 1363, 1387, 2674, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_5(buffer, 4720, 0, 3, 3712, 3748, 2419, 4219, 1387, 1411, 2728, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_8(buffer, 4828, 0, 3, 3748, 3829, 2470, 4318, 1411, 1435, 2782, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_8(buffer, 4936, 0, 3, 3829, 3910, 2521, 4417, 1435, 1459, 2836, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_9(buffer, 5044, 0, 3, 3910, 3991, 2572, 4516, 1459, 1483, 2890, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5152, 3, 1531, 1534, 2947, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5155, 3, 1534, 1537, 2950, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5158, 3, 1537, 1540, 2953, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5161, 3, 1540, 1543, 2956, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5164, 3, 1543, 1546, 2959, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5167, 3, 1546, 1549, 2962, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5170, 3, 1549, 1552, 2965, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5173, 3, 1552, 1555, 2968, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5176, 3, 1555, 1558, 2971, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5179, 0, 2944, 5152, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5182, 0, 2947, 5155, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5185, 0, 2950, 5158, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5188, 0, 2953, 5161, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5191, 0, 2956, 5164, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5194, 0, 2959, 5167, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5197, 0, 2962, 5170, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5200, 0, 2965, 5173, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5203, 0, 2968, 5176, ncols, p);

            compute_prim_dg_electron_repulsion_27(buffer, 5206, 3, 2974, 1597, 1600, 3007, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5209, 3, 2977, 1600, 1603, 3013, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5218, 3, 2980, 1603, 1606, 3019, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5227, 3, 2983, 1606, 1609, 3025, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5236, 3, 2986, 1609, 1612, 3031, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5245, 3, 2989, 1612, 1615, 3037, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5254, 3, 2992, 1615, 1618, 3043, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5263, 3, 2995, 1618, 1621, 3049, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 5272, 3, 2998, 1621, 1624, 3055, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_38(buffer, 5275, 0, 3, 3007, 5209, 1630, 1639, 3067, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_55(buffer, 5290, 0, 3, 3013, 5218, 1639, 1648, 3079, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 5305, 0, 3, 3019, 5227, 1648, 1657, 3094, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 5326, 0, 3, 3025, 5236, 1657, 1666, 3109, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 5347, 0, 3, 3031, 5245, 1666, 1675, 3124, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 5368, 0, 3, 3037, 5254, 1675, 1684, 3139, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 5389, 0, 3, 3043, 5263, 1684, 1693, 3154, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_54(buffer, 5410, 0, 3, 3049, 5272, 1693, 1702, 3169, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_62(buffer, 5419, 0, 3, 5206, 5209, 3067, 5290, 1720, 1735, 3214, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_63(buffer, 5440, 0, 3, 5209, 5218, 3079, 5305, 1735, 1750, 3238, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_64(buffer, 5473, 0, 3, 5218, 5227, 3094, 5326, 1750, 1765, 3262, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 5506, 0, 3, 5227, 5236, 3109, 5347, 1765, 1786, 3292, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 5545, 0, 3, 5236, 5245, 3124, 5368, 1786, 1807, 3322, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_66(buffer, 5584, 0, 3, 5245, 5254, 3139, 5389, 1807, 1828, 3352, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_67(buffer, 5617, 0, 3, 5254, 5263, 3154, 5410, 1828, 1849, 3376, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_53(buffer, 5638, 0, 3, 5275, 5290, 3214, 5440, 1879, 1900, 3421, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_54(buffer, 5678, 0, 3, 5290, 5305, 3238, 5473, 1900, 1921, 3451, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_55(buffer, 5718, 0, 3, 5305, 5326, 3262, 5506, 1921, 1954, 3502, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_56(buffer, 5791, 0, 3, 5326, 5347, 3292, 5545, 1954, 1987, 3553, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_57(buffer, 5861, 0, 3, 5347, 5368, 3322, 5584, 1987, 2020, 3604, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_58(buffer, 5925, 0, 3, 5368, 5389, 3352, 5617, 2020, 2053, 3649, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_34(buffer, 5959, 0, 3, 5419, 5440, 3421, 5678, 2101, 2128, 3748, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_35(buffer, 6007, 0, 3, 5440, 5473, 3451, 5718, 2128, 2170, 3829, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_36(buffer, 6130, 0, 3, 5473, 5506, 3502, 5791, 2170, 2212, 3910, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_37(buffer, 6259, 0, 3, 5506, 5545, 3553, 5861, 2212, 2254, 3991, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_38(buffer, 6382, 0, 3, 5545, 5584, 3604, 5925, 2254, 2296, 4072, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_16(buffer, 6454, 0, 3, 5638, 5678, 3748, 6007, 2371, 2419, 4219, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_17(buffer, 6604, 0, 3, 5678, 5718, 3829, 6130, 2419, 2470, 4318, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_18(buffer, 6757, 0, 3, 5718, 5791, 3910, 6259, 2470, 2521, 4417, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_19(buffer, 6913, 0, 3, 5791, 5861, 3991, 6382, 2521, 2572, 4516, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_5(buffer, 7060, 0, 3, 5959, 6007, 4219, 6604, 2674, 2728, 4828, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_6(buffer, 7228, 0, 3, 6007, 6130, 4318, 6757, 2728, 2782, 4936, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_7(buffer, 7396, 0, 3, 6130, 6259, 4417, 6913, 2782, 2836, 5044, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7564, 3, 2944, 2947, 5155, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7567, 3, 2947, 2950, 5158, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7570, 3, 2950, 2953, 5161, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7573, 3, 2953, 2956, 5164, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7576, 3, 2956, 2959, 5167, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7579, 3, 2959, 2962, 5170, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7582, 3, 2962, 2965, 5173, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 7585, 3, 2965, 2968, 5176, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7588, 0, 5155, 7567, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7591, 0, 5158, 7570, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7594, 0, 5161, 7573, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7597, 0, 5164, 7576, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7600, 0, 5167, 7579, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7603, 0, 5170, 7582, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 7606, 0, 5173, 7585, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 7609, 3, 5179, 3001, 3004, 5206, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_32(buffer, 7612, 3, 5182, 3004, 3007, 5209, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7615, 3, 5185, 3007, 3013, 5218, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7627, 3, 5188, 3013, 3019, 5227, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7639, 3, 5191, 3019, 3025, 5236, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7651, 3, 5194, 3025, 3031, 5245, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7663, 3, 5197, 3031, 3037, 5254, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 7675, 3, 5200, 3037, 3043, 5263, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 7687, 3, 5203, 3043, 3049, 5272, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_34(buffer, 7690, 0, 3, 5209, 7615, 3058, 3067, 5290, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_59(buffer, 7699, 0, 3, 5218, 7627, 3067, 3079, 5305, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 7717, 0, 3, 5227, 7639, 3079, 3094, 5326, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 7744, 0, 3, 5236, 7651, 3094, 3109, 5347, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 7771, 0, 3, 5245, 7663, 3109, 3124, 5368, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 7798, 0, 3, 5254, 7675, 3124, 3139, 5389, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_58(buffer, 7825, 0, 3, 5263, 7687, 3139, 3154, 5410, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_24(buffer, 7834, 0, 3, 7609, 7612, 5275, 7690, 3178, 3196, 5419, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_57(buffer, 7858, 0, 3, 7612, 7615, 5290, 7699, 3196, 3214, 5440, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_58(buffer, 7882, 0, 3, 7615, 7627, 5305, 7717, 3214, 3238, 5473, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_59(buffer, 7924, 0, 3, 7627, 7639, 5326, 7744, 3238, 3262, 5506, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_60(buffer, 7966, 0, 3, 7639, 7651, 5347, 7771, 3262, 3292, 5545, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_61(buffer, 8014, 0, 3, 7651, 7663, 5368, 7798, 3292, 3322, 5584, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_62(buffer, 8056, 0, 3, 7663, 7675, 5389, 7825, 3322, 3352, 5617, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_40(buffer, 8080, 0, 3, 7690, 7699, 5440, 7882, 3394, 3421, 5678, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_41(buffer, 8121, 0, 3, 7699, 7717, 5473, 7924, 3421, 3451, 5718, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_42(buffer, 8174, 0, 3, 7717, 7744, 5506, 7966, 3451, 3502, 5791, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_43(buffer, 8272, 0, 3, 7744, 7771, 5545, 8014, 3502, 3553, 5861, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_44(buffer, 8358, 0, 3, 7771, 7798, 5584, 8056, 3553, 3604, 5925, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_24(buffer, 8402, 0, 3, 7834, 7858, 5638, 8080, 3676, 3712, 5959, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_25(buffer, 8462, 0, 3, 7858, 7882, 5678, 8121, 3712, 3748, 6007, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_26(buffer, 8522, 0, 3, 7882, 7924, 5718, 8174, 3748, 3829, 6130, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_27(buffer, 8744, 0, 3, 7924, 7966, 5791, 8272, 3829, 3910, 6259, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_28(buffer, 8933, 0, 3, 7966, 8014, 5861, 8358, 3910, 3991, 6382, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_11(buffer, 9038, 0, 3, 8080, 8121, 6007, 8522, 4123, 4219, 6604, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_12(buffer, 9251, 0, 3, 8121, 8174, 6130, 8744, 4219, 4318, 6757, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_13(buffer, 9534, 0, 3, 8174, 8272, 6259, 8933, 4318, 4417, 6913, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_2(buffer, 9771, 0, 3, 8402, 8462, 6454, 9038, 4612, 4720, 7060, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_3(buffer, 10023, 0, 3, 8462, 8522, 6604, 9251, 4720, 4828, 7228, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_4(buffer, 10275, 0, 3, 8522, 8744, 6757, 9534, 4828, 4936, 7396, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10587, 3, 5152, 5155, 7567, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10590, 3, 5155, 5158, 7570, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10593, 3, 5158, 5161, 7573, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10596, 3, 5161, 5164, 7576, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10599, 3, 5164, 5167, 7579, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10602, 3, 5167, 5170, 7582, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 10605, 3, 5170, 5173, 7585, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10608, 0, 7564, 10587, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10611, 0, 7567, 10590, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10614, 0, 7570, 10593, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10617, 0, 7573, 10596, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10620, 0, 7576, 10599, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10623, 0, 7579, 10602, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 10626, 0, 7582, 10605, ncols, p);

            compute_prim_di_electron_repulsion_29(buffer, 10629, 3, 7588, 5206, 5209, 7615, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 10632, 3, 7591, 5209, 5218, 7627, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 10644, 3, 7594, 5218, 5227, 7639, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 10656, 3, 7597, 5227, 5236, 7651, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 10668, 3, 7600, 5236, 5245, 7663, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_30(buffer, 10680, 3, 7603, 5245, 5254, 7675, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 10692, 3, 7606, 5254, 5263, 7687, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_23(buffer, 10695, 0, 3, 7615, 10632, 5275, 5290, 7699, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_46(buffer, 10716, 0, 3, 7627, 10644, 5290, 5305, 7717, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_36(buffer, 10737, 0, 3, 7639, 10656, 5305, 5326, 7744, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_36(buffer, 10767, 0, 3, 7651, 10668, 5326, 5347, 7771, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_36(buffer, 10797, 0, 3, 7663, 10680, 5347, 5368, 7798, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_47(buffer, 10827, 0, 3, 7675, 10692, 5368, 5389, 7825, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_38(buffer, 10836, 0, 3, 10629, 10632, 7699, 10716, 5419, 5440, 7882, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_39(buffer, 10863, 0, 3, 10632, 10644, 7717, 10737, 5440, 5473, 7924, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_40(buffer, 10911, 0, 3, 10644, 10656, 7744, 10767, 5473, 5506, 7966, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_41(buffer, 10959, 0, 3, 10656, 10668, 7771, 10797, 5506, 5545, 8014, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_42(buffer, 11007, 0, 3, 10668, 10680, 7798, 10827, 5545, 5584, 8056, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_25(buffer, 11034, 0, 3, 10695, 10716, 7882, 10863, 5638, 5678, 8121, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_26(buffer, 11103, 0, 3, 10716, 10737, 7924, 10911, 5678, 5718, 8174, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_27(buffer, 11172, 0, 3, 10737, 10767, 7966, 10959, 5718, 5791, 8272, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_28(buffer, 11274, 0, 3, 10767, 10797, 8014, 11007, 5791, 5861, 8358, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_14(buffer, 11331, 0, 3, 10836, 10863, 8121, 11103, 5959, 6007, 8522, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_15(buffer, 11430, 0, 3, 10863, 10911, 8174, 11172, 6007, 6130, 8744, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_16(buffer, 11696, 0, 3, 10911, 10959, 8272, 11274, 6130, 6259, 8933, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_6(buffer, 11846, 0, 3, 11034, 11103, 8522, 11430, 6454, 6604, 9251, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_7(buffer, 12486, 0, 3, 11103, 11172, 8744, 11696, 6604, 6757, 9534, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_1(buffer, 12884, 0, 3, 11331, 11430, 9251, 12486, 7060, 7228, 10275, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_5(buffer, 13664, 3, 10608, 7609, 7612, 10629, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_16(buffer, 13667, 3, 10611, 7612, 7615, 10632, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 13670, 3, 10614, 7615, 7627, 10644, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 13673, 3, 10617, 7627, 7639, 10656, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 13676, 3, 10620, 7639, 7651, 10668, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_17(buffer, 13679, 3, 10623, 7651, 7663, 10680, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_24(buffer, 13682, 3, 10626, 7663, 7675, 10692, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_11(buffer, 13685, 0, 3, 10632, 13670, 7690, 7699, 10716, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_27(buffer, 13694, 0, 3, 10644, 13673, 7699, 7717, 10737, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_19(buffer, 13703, 0, 3, 10656, 13676, 7717, 7744, 10767, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_19(buffer, 13712, 0, 3, 10668, 13679, 7744, 7771, 10797, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_28(buffer, 13721, 0, 3, 10680, 13682, 7771, 7798, 10827, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_5(buffer, 13730, 0, 3, 13664, 13667, 10695, 13685, 7834, 7858, 10836, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_18(buffer, 13757, 0, 3, 13667, 13670, 10716, 13694, 7858, 7882, 10863, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_19(buffer, 13784, 0, 3, 13670, 13673, 10737, 13703, 7882, 7924, 10911, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_20(buffer, 13811, 0, 3, 13673, 13676, 10767, 13712, 7924, 7966, 10959, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_21(buffer, 13838, 0, 3, 13676, 13679, 10797, 13721, 7966, 8014, 11007, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_11(buffer, 13865, 0, 3, 13685, 13694, 10863, 13784, 8080, 8121, 11103, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_12(buffer, 13925, 0, 3, 13694, 13703, 10911, 13811, 8121, 8174, 11172, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_13(buffer, 13985, 0, 3, 13703, 13712, 10959, 13838, 8174, 8272, 11274, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_5(buffer, 14045, 0, 3, 13730, 13757, 11034, 13865, 8402, 8462, 11331, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_6(buffer, 14153, 0, 3, 13757, 13784, 11103, 13925, 8462, 8522, 11430, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_7(buffer, 14261, 0, 3, 13784, 13811, 11172, 13985, 8522, 8744, 11696, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_2(buffer, 14420, 0, 3, 13865, 13925, 11430, 14261, 9038, 9251, 12486, ncols, alpha, beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 14879, 0, 3, 14045, 14153, 11846, 14420, 9771, 10023, 12884, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 16499, 14879, 1620, ncols);
        }
    }

    simdtrf::transform_lk(values, nvalues, buffer, 16499, nmax);
}

}  // namespace simdt2ceri
