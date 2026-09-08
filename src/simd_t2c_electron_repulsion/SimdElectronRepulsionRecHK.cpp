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


#include "SimdElectronRepulsionRecHK.hpp"

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
#include "SimdTransformHK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(6439, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 64, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 67, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 70, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 73, 0, 14, 15, 43, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 76, 0, 15, 16, 46, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 79, 0, 16, 17, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 19, 22, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 88, 0, 22, 25, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 94, 0, 25, 28, 61, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 28, 31, 64, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 106, 0, 31, 34, 67, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 112, 0, 34, 37, 70, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 118, 0, 37, 40, 73, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 124, 0, 40, 43, 76, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 130, 0, 43, 46, 79, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 136, 0, 52, 55, 88, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 142, 0, 55, 58, 94, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 151, 0, 58, 61, 100, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 160, 0, 61, 64, 106, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 169, 0, 64, 67, 112, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 178, 0, 67, 70, 118, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 187, 0, 70, 73, 124, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 196, 0, 73, 76, 130, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 205, 0, 82, 88, 142, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 214, 0, 88, 94, 151, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 223, 0, 94, 100, 160, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 232, 0, 100, 106, 169, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 241, 0, 106, 112, 178, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 250, 0, 112, 118, 187, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 259, 0, 118, 124, 196, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 268, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 271, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 274, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 277, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 280, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 283, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 286, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 289, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 292, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 295, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 298, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 301, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 304, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 307, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 310, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 313, 3, 16, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 316, 3, 22, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 319, 3, 25, 58, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 322, 3, 28, 61, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 331, 3, 31, 64, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 340, 3, 34, 67, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 349, 3, 37, 70, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 358, 3, 40, 73, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 367, 3, 43, 76, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 376, 3, 46, 79, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 379, 3, 52, 82, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 382, 3, 55, 88, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 385, 3, 58, 94, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 394, 3, 61, 100, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 403, 3, 64, 106, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 412, 3, 67, 112, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 421, 3, 70, 118, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 430, 3, 73, 124, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 439, 3, 76, 130, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 448, 3, 88, 142, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 460, 3, 94, 151, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 472, 3, 100, 160, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 484, 3, 106, 169, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 496, 3, 112, 178, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 508, 3, 118, 187, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 520, 3, 124, 196, ncols, p);

            compute_prim_hp_electron_repulsion_2(buffer, 532, 3, 136, 205, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 547, 3, 142, 214, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 562, 3, 151, 223, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 577, 3, 160, 232, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 592, 3, 169, 241, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 607, 3, 178, 250, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 622, 3, 187, 259, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 637, 3, 9, 10, 271, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 640, 3, 10, 11, 274, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 643, 3, 11, 12, 277, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 646, 3, 12, 13, 280, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 649, 3, 13, 14, 283, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 652, 3, 14, 15, 286, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 655, 3, 15, 16, 289, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 658, 0, 268, 637, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 661, 0, 271, 640, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 664, 0, 274, 643, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 667, 0, 277, 646, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 670, 0, 280, 649, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 673, 0, 283, 652, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 676, 0, 286, 655, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 679, 3, 292, 52, 55, 319, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 682, 3, 295, 55, 58, 322, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 685, 3, 298, 58, 61, 331, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 694, 3, 301, 61, 64, 340, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 703, 3, 304, 64, 67, 349, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 712, 3, 307, 67, 70, 358, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 721, 3, 310, 70, 73, 367, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 730, 3, 313, 73, 76, 376, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 733, 0, 3, 319, 682, 82, 88, 385, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 742, 0, 3, 322, 685, 88, 94, 394, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 757, 0, 3, 331, 694, 94, 100, 403, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 772, 0, 3, 340, 703, 100, 106, 412, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 787, 0, 3, 349, 712, 106, 112, 421, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 802, 0, 3, 358, 721, 112, 118, 430, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 817, 0, 3, 367, 730, 118, 124, 439, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_14(buffer, 832, 0, 3, 679, 682, 385, 742, 136, 142, 460, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 853, 0, 3, 682, 685, 394, 757, 142, 151, 472, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 877, 0, 3, 685, 694, 403, 772, 151, 160, 484, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 901, 0, 3, 694, 703, 412, 787, 160, 169, 496, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 925, 0, 3, 703, 712, 421, 802, 169, 178, 508, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 949, 0, 3, 712, 721, 430, 817, 178, 187, 520, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_5(buffer, 973, 0, 3, 733, 742, 460, 853, 205, 214, 562, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1000, 0, 3, 742, 757, 472, 877, 214, 223, 577, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1027, 0, 3, 757, 772, 484, 901, 223, 232, 592, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1054, 0, 3, 772, 787, 496, 925, 232, 241, 607, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1081, 0, 3, 787, 802, 508, 949, 241, 250, 622, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1108, 3, 268, 271, 640, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1111, 3, 271, 274, 643, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1114, 3, 274, 277, 646, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1117, 3, 277, 280, 649, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1120, 3, 280, 283, 652, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1123, 3, 283, 286, 655, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1126, 0, 637, 1108, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1129, 0, 640, 1111, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1132, 0, 643, 1114, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1135, 0, 646, 1117, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1138, 0, 649, 1120, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1141, 0, 652, 1123, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1144, 3, 658, 316, 319, 682, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 1147, 3, 661, 319, 322, 685, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1150, 3, 664, 322, 331, 694, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1168, 3, 667, 331, 340, 703, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1186, 3, 670, 340, 349, 712, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1204, 3, 673, 349, 358, 721, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 1222, 3, 676, 358, 367, 730, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1225, 0, 3, 679, 1144, 379, 382, 733, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_21(buffer, 1234, 0, 3, 682, 1147, 382, 385, 742, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1243, 0, 3, 685, 1150, 385, 394, 757, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1270, 0, 3, 694, 1168, 394, 403, 772, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1297, 0, 3, 703, 1186, 403, 412, 787, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1324, 0, 3, 712, 1204, 412, 421, 802, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_20(buffer, 1351, 0, 3, 721, 1222, 421, 430, 817, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_20(buffer, 1375, 0, 3, 1144, 1147, 742, 1243, 448, 460, 853, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_21(buffer, 1417, 0, 3, 1147, 1150, 757, 1270, 460, 472, 877, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1462, 0, 3, 1150, 1168, 772, 1297, 472, 484, 901, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1507, 0, 3, 1168, 1186, 787, 1324, 484, 496, 925, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_19(buffer, 1552, 0, 3, 1186, 1204, 802, 1351, 496, 508, 949, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_2(buffer, 1594, 0, 3, 1225, 1234, 832, 1375, 532, 547, 973, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_5(buffer, 1648, 0, 3, 1234, 1243, 853, 1417, 547, 562, 1000, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_8(buffer, 1702, 0, 3, 1243, 1270, 877, 1462, 562, 577, 1027, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_8(buffer, 1756, 0, 3, 1270, 1297, 901, 1507, 577, 592, 1054, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_9(buffer, 1810, 0, 3, 1297, 1324, 925, 1552, 592, 607, 1081, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1864, 3, 637, 640, 1111, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1867, 3, 640, 643, 1114, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1870, 3, 643, 646, 1117, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1873, 3, 646, 649, 1120, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1876, 3, 649, 652, 1123, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1879, 0, 1108, 1864, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1882, 0, 1111, 1867, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1885, 0, 1114, 1870, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1888, 0, 1117, 1873, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1891, 0, 1120, 1876, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1894, 3, 1126, 679, 682, 1147, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_14(buffer, 1897, 3, 1129, 682, 685, 1150, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_17(buffer, 1900, 3, 1132, 685, 694, 1168, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1930, 3, 1135, 694, 703, 1186, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1957, 3, 1138, 703, 712, 1204, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_22(buffer, 1984, 3, 1141, 712, 721, 1222, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_22(buffer, 1987, 0, 3, 1147, 1897, 733, 742, 1243, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_11(buffer, 1996, 0, 3, 1150, 1900, 742, 757, 1270, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_26(buffer, 2035, 0, 3, 1168, 1930, 757, 772, 1297, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 2080, 0, 3, 1186, 1957, 772, 787, 1324, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_25(buffer, 2119, 0, 3, 1204, 1984, 787, 802, 1351, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_16(buffer, 2152, 0, 3, 1894, 1897, 1243, 1996, 832, 853, 1417, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_17(buffer, 2215, 0, 3, 1897, 1900, 1270, 2035, 853, 877, 1462, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_18(buffer, 2281, 0, 3, 1900, 1930, 1297, 2080, 877, 901, 1507, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_19(buffer, 2350, 0, 3, 1930, 1957, 1324, 2119, 901, 925, 1552, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_5(buffer, 2410, 0, 3, 1987, 1996, 1417, 2215, 973, 1000, 1702, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_6(buffer, 2491, 0, 3, 1996, 2035, 1462, 2281, 1000, 1027, 1756, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_7(buffer, 2572, 0, 3, 2035, 2080, 1507, 2350, 1027, 1054, 1810, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 2653, 3, 1108, 1111, 1867, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 2656, 3, 1111, 1114, 1870, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 2659, 3, 1114, 1117, 1873, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 2662, 3, 1117, 1120, 1876, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2665, 0, 1864, 2653, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2668, 0, 1867, 2656, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2671, 0, 1870, 2659, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2674, 0, 1873, 2662, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 2677, 3, 1879, 1144, 1147, 1897, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_24(buffer, 2680, 3, 1882, 1147, 1150, 1900, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_25(buffer, 2683, 0, 3, 1885, 2671, 1150, 1168, 1930, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_26(buffer, 2722, 3, 1888, 1168, 1186, 1957, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_23(buffer, 2758, 3, 1891, 1186, 1204, 1984, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_14(buffer, 2761, 0, 3, 1894, 2677, 1225, 1234, 1987, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_20(buffer, 2770, 0, 3, 1897, 2680, 1234, 1243, 1996, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_21(buffer, 2779, 0, 3, 1900, 2683, 1243, 1270, 2035, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_22(buffer, 2857, 0, 3, 1930, 2722, 1270, 1297, 2080, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_23(buffer, 2920, 0, 3, 1957, 2758, 1297, 1324, 2119, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_11(buffer, 2965, 0, 3, 2677, 2680, 1996, 2779, 1375, 1417, 2215, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_12(buffer, 3049, 0, 3, 2680, 2683, 2035, 2857, 1417, 1462, 2281, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_13(buffer, 3170, 0, 3, 2683, 2722, 2080, 2920, 1462, 1507, 2350, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_2(buffer, 3260, 0, 3, 2761, 2770, 2152, 2965, 1594, 1648, 2410, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_3(buffer, 3377, 0, 3, 2770, 2779, 2215, 3049, 1648, 1702, 2491, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_4(buffer, 3494, 0, 3, 2779, 2857, 2281, 3170, 1702, 1756, 2572, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 3635, 3, 1864, 1867, 2656, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 3638, 3, 1867, 1870, 2659, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 3641, 3, 1870, 1873, 2662, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 3644, 0, 2653, 3635, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 3647, 0, 2656, 3638, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 3650, 0, 2659, 3641, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 3653, 3, 2665, 1894, 1897, 2680, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_20(buffer, 3656, 3, 2668, 1897, 1900, 2683, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_21(buffer, 3659, 3, 2671, 1900, 1930, 2722, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_22(buffer, 3695, 3, 2674, 1930, 1957, 2758, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_13(buffer, 3698, 0, 3, 2680, 3656, 1987, 1996, 2779, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_14(buffer, 3707, 0, 3, 2683, 3659, 1996, 2035, 2857, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_15(buffer, 3811, 0, 3, 2722, 3695, 2035, 2080, 2920, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_6(buffer, 3871, 0, 3, 3653, 3656, 2779, 3707, 2152, 2215, 3049, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_7(buffer, 4136, 0, 3, 3656, 3659, 2857, 3811, 2215, 2281, 3170, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_1(buffer, 4294, 0, 3, 3698, 3707, 3049, 4136, 2410, 2491, 3494, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_5(buffer, 4651, 3, 3644, 2677, 2680, 3656, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_11(buffer, 4654, 3, 3647, 2680, 2683, 3659, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_12(buffer, 4657, 3, 3650, 2683, 2722, 3695, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_5(buffer, 4660, 0, 3, 3653, 4651, 2761, 2770, 3698, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_6(buffer, 4669, 0, 3, 3656, 4654, 2770, 2779, 3707, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_7(buffer, 4678, 0, 3, 3659, 4657, 2779, 2857, 3811, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_2(buffer, 4738, 0, 3, 4651, 4654, 3707, 4678, 2965, 3049, 4136, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 4927, 0, 3, 4660, 4669, 3871, 4738, 3260, 3377, 4294, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5683, 4927, 756, ncols);
        }
    }

    simdtrf::transform_hk(values, nvalues, buffer, 5683, nmax);
}

}  // namespace simdt2ceri
