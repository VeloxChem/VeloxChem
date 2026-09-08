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


#include "SimdElectronRepulsionRecIL.hpp"

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
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFL.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGL.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHL.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIL.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformIL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_il_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_il_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(13524, nvalues);

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

            compute_prim_is_electron_repulsion_1(buffer, 373, 0, 166, 175, 265, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 385, 0, 175, 184, 277, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 397, 0, 184, 193, 289, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 409, 0, 193, 202, 301, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 421, 0, 202, 211, 313, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 433, 0, 211, 220, 325, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 445, 0, 220, 229, 337, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 457, 0, 229, 238, 349, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 469, 0, 238, 247, 361, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 481, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 484, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 487, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 490, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 493, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 496, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 499, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 502, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 505, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 508, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 511, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 514, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 517, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 520, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 523, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 526, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 529, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 532, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 535, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 538, 3, 17, 49, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 541, 3, 18, 52, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 544, 3, 19, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 547, 3, 25, 67, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 550, 3, 28, 70, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 553, 3, 31, 73, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 556, 3, 34, 76, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 565, 3, 37, 79, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 574, 3, 40, 82, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 583, 3, 43, 85, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 592, 3, 46, 88, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 601, 3, 49, 91, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 610, 3, 52, 94, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 613, 3, 55, 97, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 616, 3, 67, 106, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 619, 3, 70, 112, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 622, 3, 73, 118, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 631, 3, 76, 124, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 640, 3, 79, 130, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 649, 3, 82, 136, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 658, 3, 85, 142, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 667, 3, 88, 148, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 676, 3, 91, 154, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 685, 3, 94, 160, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 688, 3, 106, 184, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 691, 3, 112, 193, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 703, 3, 118, 202, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 715, 3, 124, 211, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 727, 3, 130, 220, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 739, 3, 136, 229, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 751, 3, 142, 238, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 763, 3, 148, 247, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 775, 3, 154, 256, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 784, 3, 184, 277, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 799, 3, 193, 289, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 814, 3, 202, 301, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 829, 3, 211, 313, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 844, 3, 220, 325, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 859, 3, 229, 337, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 874, 3, 238, 349, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 889, 3, 247, 361, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 904, 3, 277, 397, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 922, 3, 289, 409, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 940, 3, 301, 421, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 958, 3, 313, 433, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 976, 3, 325, 445, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 994, 3, 337, 457, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 1012, 3, 349, 469, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1030, 3, 9, 10, 484, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1033, 3, 10, 11, 487, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1036, 3, 11, 12, 490, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1039, 3, 12, 13, 493, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1042, 3, 13, 14, 496, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1045, 3, 14, 15, 499, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1048, 3, 15, 16, 502, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1051, 3, 16, 17, 505, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1054, 3, 17, 18, 508, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1057, 3, 18, 19, 511, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1060, 0, 484, 1033, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1063, 0, 487, 1036, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1066, 0, 490, 1039, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1069, 0, 493, 1042, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1072, 0, 496, 1045, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1075, 0, 499, 1048, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1078, 0, 502, 1051, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1081, 0, 505, 1054, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1084, 0, 508, 1057, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1087, 3, 514, 61, 64, 547, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1090, 3, 517, 64, 67, 550, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1093, 3, 520, 67, 70, 553, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 1096, 3, 523, 70, 73, 556, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 1099, 3, 526, 73, 76, 565, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 1108, 3, 529, 76, 79, 574, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 1117, 3, 532, 79, 82, 583, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 1126, 3, 535, 82, 85, 592, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 1135, 3, 538, 85, 88, 601, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1144, 3, 541, 88, 91, 610, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1147, 3, 544, 91, 94, 613, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1150, 0, 3, 550, 1093, 100, 106, 619, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 1159, 0, 3, 553, 1096, 106, 112, 622, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1168, 0, 3, 556, 1099, 112, 118, 631, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1183, 0, 3, 565, 1108, 118, 124, 640, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1198, 0, 3, 574, 1117, 124, 130, 649, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1213, 0, 3, 583, 1126, 130, 136, 658, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 1228, 0, 3, 592, 1135, 136, 142, 667, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 1243, 0, 3, 601, 1144, 142, 148, 676, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1258, 0, 3, 610, 1147, 148, 154, 685, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1267, 0, 3, 1087, 1090, 616, 1150, 166, 175, 688, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1282, 0, 3, 1090, 1093, 619, 1159, 175, 184, 691, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1297, 0, 3, 1093, 1096, 622, 1168, 184, 193, 703, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 1321, 0, 3, 1096, 1099, 631, 1183, 193, 202, 715, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1345, 0, 3, 1099, 1108, 640, 1198, 202, 211, 727, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1369, 0, 3, 1108, 1117, 649, 1213, 211, 220, 739, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1393, 0, 3, 1117, 1126, 658, 1228, 220, 229, 751, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1417, 0, 3, 1126, 1135, 667, 1243, 229, 238, 763, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_29(buffer, 1441, 0, 3, 1135, 1144, 676, 1258, 238, 247, 775, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1462, 0, 3, 1150, 1159, 691, 1297, 265, 277, 799, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1495, 0, 3, 1159, 1168, 703, 1321, 277, 289, 814, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1528, 0, 3, 1168, 1183, 715, 1345, 289, 301, 829, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1561, 0, 3, 1183, 1198, 727, 1369, 301, 313, 844, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1594, 0, 3, 1198, 1213, 739, 1393, 313, 325, 859, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1627, 0, 3, 1213, 1228, 751, 1417, 325, 337, 874, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1660, 0, 3, 1228, 1243, 763, 1441, 337, 349, 889, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_2(buffer, 1693, 0, 3, 1267, 1282, 784, 1462, 373, 385, 904, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_5(buffer, 1729, 0, 3, 1282, 1297, 799, 1495, 385, 397, 922, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1765, 0, 3, 1297, 1321, 814, 1528, 397, 409, 940, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1801, 0, 3, 1321, 1345, 829, 1561, 409, 421, 958, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1837, 0, 3, 1345, 1369, 844, 1594, 421, 433, 976, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1873, 0, 3, 1369, 1393, 859, 1627, 433, 445, 994, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1909, 0, 3, 1393, 1417, 874, 1660, 445, 457, 1012, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1945, 3, 481, 484, 1033, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1948, 3, 484, 487, 1036, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1951, 3, 487, 490, 1039, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1954, 3, 490, 493, 1042, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1957, 3, 493, 496, 1045, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1960, 3, 496, 499, 1048, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1963, 3, 499, 502, 1051, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1966, 3, 502, 505, 1054, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1969, 3, 505, 508, 1057, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1972, 0, 1030, 1945, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1975, 0, 1033, 1948, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1978, 0, 1036, 1951, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1981, 0, 1039, 1954, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1984, 0, 1042, 1957, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1987, 0, 1045, 1960, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1990, 0, 1048, 1963, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1993, 0, 1051, 1966, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1996, 0, 1054, 1969, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1999, 3, 1060, 547, 550, 1093, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2002, 3, 1063, 550, 553, 1096, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_25(buffer, 2008, 3, 1066, 553, 556, 1099, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 2014, 3, 1069, 556, 565, 1108, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 2026, 3, 1072, 565, 574, 1117, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 2038, 3, 1075, 574, 583, 1126, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 2050, 3, 1078, 583, 592, 1135, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_27(buffer, 2062, 3, 1081, 592, 601, 1144, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_23(buffer, 2068, 3, 1084, 601, 610, 1147, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 2071, 0, 3, 1093, 2002, 616, 619, 1159, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_32(buffer, 2083, 0, 3, 1096, 2008, 619, 622, 1168, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 2095, 0, 3, 1099, 2014, 622, 631, 1183, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 2122, 0, 3, 1108, 2026, 631, 640, 1198, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 2149, 0, 3, 1117, 2038, 640, 649, 1213, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 2176, 0, 3, 1126, 2050, 649, 658, 1228, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_34(buffer, 2203, 0, 3, 1135, 2062, 658, 667, 1243, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_30(buffer, 2230, 0, 3, 1144, 2068, 667, 676, 1258, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_30(buffer, 2239, 0, 3, 1999, 2002, 1159, 2083, 688, 691, 1297, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_37(buffer, 2257, 0, 3, 2002, 2008, 1168, 2095, 691, 703, 1321, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_38(buffer, 2302, 0, 3, 2008, 2014, 1183, 2122, 703, 715, 1345, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 2347, 0, 3, 2014, 2026, 1198, 2149, 715, 727, 1369, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 2392, 0, 3, 2026, 2038, 1213, 2176, 727, 739, 1393, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 2437, 0, 3, 2038, 2050, 1228, 2203, 739, 751, 1417, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_33(buffer, 2482, 0, 3, 2050, 2062, 1243, 2230, 751, 763, 1441, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_16(buffer, 2515, 0, 3, 2071, 2083, 1297, 2257, 784, 799, 1495, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_21(buffer, 2575, 0, 3, 2083, 2095, 1321, 2302, 799, 814, 1528, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_22(buffer, 2638, 0, 3, 2095, 2122, 1345, 2347, 814, 829, 1561, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_22(buffer, 2701, 0, 3, 2122, 2149, 1369, 2392, 829, 844, 1594, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_22(buffer, 2764, 0, 3, 2149, 2176, 1393, 2437, 844, 859, 1627, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_19(buffer, 2827, 0, 3, 2176, 2203, 1417, 2482, 859, 874, 1660, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_5(buffer, 2887, 0, 3, 2239, 2257, 1495, 2575, 904, 922, 1765, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_8(buffer, 2959, 0, 3, 2257, 2302, 1528, 2638, 922, 940, 1801, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_8(buffer, 3031, 0, 3, 2302, 2347, 1561, 2701, 940, 958, 1837, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_8(buffer, 3103, 0, 3, 2347, 2392, 1594, 2764, 958, 976, 1873, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_9(buffer, 3175, 0, 3, 2392, 2437, 1627, 2827, 976, 994, 1909, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3247, 3, 1030, 1033, 1948, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3250, 3, 1033, 1036, 1951, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3253, 3, 1036, 1039, 1954, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3256, 3, 1039, 1042, 1957, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3259, 3, 1042, 1045, 1960, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3262, 3, 1045, 1048, 1963, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3265, 3, 1048, 1051, 1966, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 3268, 3, 1051, 1054, 1969, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3271, 0, 1948, 3250, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3274, 0, 1951, 3253, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3277, 0, 1954, 3256, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3280, 0, 1957, 3259, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3283, 0, 1960, 3262, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3286, 0, 1963, 3265, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 3289, 0, 1966, 3268, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 3292, 3, 1972, 1087, 1090, 1999, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 3295, 3, 1975, 1090, 1093, 2002, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 3298, 3, 1978, 1093, 1096, 2008, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_32(buffer, 3307, 3, 1981, 1096, 1099, 2014, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 3316, 3, 1984, 1099, 1108, 2026, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 3331, 3, 1987, 1108, 1117, 2038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 3346, 3, 1990, 1117, 1126, 2050, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_34(buffer, 3361, 3, 1993, 1126, 1135, 2062, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_26(buffer, 3370, 3, 1996, 1135, 1144, 2068, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_34(buffer, 3373, 0, 3, 2002, 3298, 1150, 1159, 2083, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_39(buffer, 3382, 0, 3, 2008, 3307, 1159, 1168, 2095, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_41(buffer, 3397, 0, 3, 2014, 3316, 1168, 1183, 2122, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_41(buffer, 3436, 0, 3, 2026, 3331, 1183, 1198, 2149, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_41(buffer, 3475, 0, 3, 2038, 3346, 1198, 1213, 2176, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_37(buffer, 3514, 0, 3, 2050, 3361, 1213, 1228, 2203, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_33(buffer, 3553, 0, 3, 2062, 3370, 1228, 1243, 2230, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_29(buffer, 3562, 0, 3, 3292, 3295, 2071, 3373, 1267, 1282, 2239, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_30(buffer, 3583, 0, 3, 3295, 3298, 2083, 3382, 1282, 1297, 2257, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_37(buffer, 3604, 0, 3, 3298, 3307, 2095, 3397, 1297, 1321, 2302, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_38(buffer, 3670, 0, 3, 3307, 3316, 2122, 3436, 1321, 1345, 2347, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_36(buffer, 3736, 0, 3, 3316, 3331, 2149, 3475, 1345, 1369, 2392, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_36(buffer, 3802, 0, 3, 3331, 3346, 2176, 3514, 1369, 1393, 2437, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_33(buffer, 3868, 0, 3, 3346, 3361, 2203, 3553, 1393, 1417, 2482, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_20(buffer, 3913, 0, 3, 3373, 3382, 2257, 3604, 1462, 1495, 2575, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_21(buffer, 4001, 0, 3, 3382, 3397, 2302, 3670, 1495, 1528, 2638, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_22(buffer, 4095, 0, 3, 3397, 3436, 2347, 3736, 1528, 1561, 2701, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_22(buffer, 4189, 0, 3, 3436, 3475, 2392, 3802, 1561, 1594, 2764, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_19(buffer, 4283, 0, 3, 3475, 3514, 2437, 3868, 1594, 1627, 2827, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_2(buffer, 4371, 0, 3, 3562, 3583, 2515, 3913, 1693, 1729, 2887, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_5(buffer, 4480, 0, 3, 3583, 3604, 2575, 4001, 1729, 1765, 2959, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_8(buffer, 4589, 0, 3, 3604, 3670, 2638, 4095, 1765, 1801, 3031, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_8(buffer, 4698, 0, 3, 3670, 3736, 2701, 4189, 1801, 1837, 3103, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_9(buffer, 4807, 0, 3, 3736, 3802, 2764, 4283, 1837, 1873, 3175, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4916, 3, 1945, 1948, 3250, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4919, 3, 1948, 1951, 3253, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4922, 3, 1951, 1954, 3256, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4925, 3, 1954, 1957, 3259, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4928, 3, 1957, 1960, 3262, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4931, 3, 1960, 1963, 3265, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 4934, 3, 1963, 1966, 3268, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4937, 0, 3247, 4916, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4940, 0, 3250, 4919, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4943, 0, 3253, 4922, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4946, 0, 3256, 4925, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4949, 0, 3259, 4928, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4952, 0, 3262, 4931, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 4955, 0, 3265, 4934, ncols, p);

            compute_prim_dh_electron_repulsion_32(buffer, 4958, 3, 3271, 1999, 2002, 3298, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 4961, 3, 3274, 2002, 2008, 3307, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_37(buffer, 4973, 3, 3277, 2008, 2014, 3316, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_38(buffer, 4985, 3, 3280, 2014, 2026, 3331, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_38(buffer, 5003, 3, 3283, 2026, 2038, 3346, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_39(buffer, 5021, 3, 3286, 2038, 2050, 3361, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_31(buffer, 5033, 3, 3289, 2050, 2062, 3370, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_38(buffer, 5036, 0, 3, 3298, 4961, 2071, 2083, 3382, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_39(buffer, 5054, 0, 3, 3307, 4973, 2083, 2095, 3397, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_40(buffer, 5072, 0, 3, 3316, 4985, 2095, 2122, 3436, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_41(buffer, 5126, 0, 3, 3331, 5003, 2122, 2149, 3475, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_37(buffer, 5177, 0, 3, 3346, 5021, 2149, 2176, 3514, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_33(buffer, 5228, 0, 3, 3361, 5033, 2176, 2203, 3553, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_25(buffer, 5237, 0, 3, 4958, 4961, 3382, 5054, 2239, 2257, 3604, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_29(buffer, 5261, 0, 3, 4961, 4973, 3397, 5072, 2257, 2302, 3670, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_30(buffer, 5352, 0, 3, 4973, 4985, 3436, 5126, 2302, 2347, 3736, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_31(buffer, 5449, 0, 3, 4985, 5003, 3475, 5177, 2347, 2392, 3802, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_28(buffer, 5540, 0, 3, 5003, 5021, 3514, 5228, 2392, 2437, 3868, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_14(buffer, 5600, 0, 3, 5036, 5054, 3604, 5261, 2515, 2575, 4001, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_15(buffer, 5732, 0, 3, 5054, 5072, 3670, 5352, 2575, 2638, 4095, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_16(buffer, 5867, 0, 3, 5072, 5126, 3736, 5449, 2638, 2701, 4189, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_17(buffer, 6005, 0, 3, 5126, 5177, 3802, 5540, 2701, 2764, 4283, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_5(buffer, 6130, 0, 3, 5237, 5261, 4001, 5732, 2887, 2959, 4589, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_6(buffer, 6289, 0, 3, 5261, 5352, 4095, 5867, 2959, 3031, 4698, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_7(buffer, 6448, 0, 3, 5352, 5449, 4189, 6005, 3031, 3103, 4807, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6607, 3, 3247, 3250, 4919, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6610, 3, 3250, 3253, 4922, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6613, 3, 3253, 3256, 4925, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6616, 3, 3256, 3259, 4928, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6619, 3, 3259, 3262, 4931, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 6622, 3, 3262, 3265, 4934, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 6625, 0, 4919, 6610, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 6628, 0, 4922, 6613, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 6631, 0, 4925, 6616, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 6634, 0, 4928, 6619, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 6637, 0, 4931, 6622, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 6640, 3, 4937, 3292, 3295, 4958, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_29(buffer, 6643, 3, 4940, 3295, 3298, 4961, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 6646, 3, 4943, 3298, 3307, 4973, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_34(buffer, 6661, 3, 4946, 3307, 3316, 4985, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_35(buffer, 6676, 3, 4949, 3316, 3331, 5003, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_36(buffer, 6697, 3, 4952, 3331, 3346, 5021, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_28(buffer, 6712, 3, 4955, 3346, 3361, 5033, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_27(buffer, 6715, 0, 3, 4961, 6646, 3373, 3382, 5054, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_28(buffer, 6724, 0, 3, 4973, 6661, 3382, 3397, 5072, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_29(buffer, 6745, 0, 3, 4985, 6676, 3397, 3436, 5126, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_30(buffer, 6814, 0, 3, 5003, 6697, 3436, 3475, 5177, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_26(buffer, 6877, 0, 3, 5021, 6712, 3475, 3514, 5228, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_17(buffer, 6886, 0, 3, 6640, 6643, 5036, 6715, 3562, 3583, 5237, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_18(buffer, 6913, 0, 3, 6643, 6646, 5054, 6724, 3583, 3604, 5261, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_19(buffer, 6940, 0, 3, 6646, 6661, 5072, 6745, 3604, 3670, 5352, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_20(buffer, 7093, 0, 3, 6661, 6676, 5126, 6814, 3670, 3736, 5449, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_21(buffer, 7222, 0, 3, 6676, 6697, 5177, 6877, 3736, 3802, 5540, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_8(buffer, 7300, 0, 3, 6715, 6724, 5261, 6940, 3913, 4001, 5732, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_9(buffer, 7471, 0, 3, 6724, 6745, 5352, 7093, 4001, 4095, 5867, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_10(buffer, 7696, 0, 3, 6745, 6814, 5449, 7222, 4095, 4189, 6005, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_2(buffer, 7876, 0, 3, 6886, 6913, 5600, 7300, 4371, 4480, 6130, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_3(buffer, 8098, 0, 3, 6913, 6940, 5732, 7471, 4480, 4589, 6289, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_4(buffer, 8320, 0, 3, 6940, 7093, 5867, 7696, 4589, 4698, 6448, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 8575, 3, 4916, 4919, 6610, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 8578, 3, 4919, 4922, 6613, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 8581, 3, 4922, 4925, 6616, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 8584, 3, 4925, 4928, 6619, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 8587, 3, 4928, 4931, 6622, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_8(buffer, 8590, 0, 6607, 8575, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 8593, 0, 6610, 8578, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 8596, 0, 6613, 8581, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 8599, 0, 6616, 8584, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 8602, 0, 6619, 8587, ncols, p);

            compute_prim_dk_electron_repulsion_20(buffer, 8605, 3, 6625, 4958, 4961, 6646, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 8608, 3, 6628, 4961, 4973, 6661, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_22(buffer, 8623, 3, 6631, 4973, 4985, 6676, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_23(buffer, 8638, 3, 6634, 4985, 5003, 6697, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_19(buffer, 8653, 3, 6637, 5003, 5021, 6712, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_14(buffer, 8656, 0, 3, 6646, 8608, 5036, 5054, 6724, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_15(buffer, 8680, 0, 3, 6661, 8623, 5054, 5072, 6745, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_16(buffer, 8704, 0, 3, 6676, 8638, 5072, 5126, 6814, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_17(buffer, 8770, 0, 3, 6697, 8653, 5126, 5177, 6877, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_8(buffer, 8779, 0, 3, 8605, 8608, 6724, 8680, 5237, 5261, 6940, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_9(buffer, 8809, 0, 3, 8608, 8623, 6745, 8704, 5261, 5352, 7093, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_10(buffer, 8993, 0, 3, 8623, 8638, 6814, 8770, 5352, 5449, 7222, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_3(buffer, 9092, 0, 3, 8656, 8680, 6940, 8809, 5600, 5732, 7471, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_4(buffer, 9569, 0, 3, 8680, 8704, 7093, 8993, 5732, 5867, 7696, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_1(buffer, 9861, 0, 3, 8779, 8809, 7471, 9569, 6130, 6289, 8320, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_2(buffer, 10464, 3, 8590, 6640, 6643, 8605, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_7(buffer, 10467, 3, 8593, 6643, 6646, 8608, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 10470, 3, 8596, 6646, 6661, 8623, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_9(buffer, 10473, 3, 8599, 6661, 6676, 8638, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_10(buffer, 10476, 3, 8602, 6676, 6697, 8653, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_5(buffer, 10479, 0, 3, 8608, 10470, 6715, 6724, 8680, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_6(buffer, 10488, 0, 3, 8623, 10473, 6724, 6745, 8704, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_7(buffer, 10497, 0, 3, 8638, 10476, 6745, 6814, 8770, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_2(buffer, 10506, 0, 3, 10464, 10467, 8656, 10479, 6886, 6913, 8779, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_3(buffer, 10536, 0, 3, 10467, 10470, 8680, 10488, 6913, 6940, 8809, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_4(buffer, 10566, 0, 3, 10470, 10473, 8704, 10497, 6940, 7093, 8993, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_1(buffer, 10665, 0, 3, 10479, 10488, 8809, 10566, 7300, 7471, 9569, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 11004, 0, 3, 10506, 10536, 9092, 10665, 7876, 8098, 9861, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 12264, 11004, 1260, ncols);
        }
    }

    simdtrf::transform_il(values, nvalues, buffer, 12264, nmax);
}

}  // namespace simdt2ceri
