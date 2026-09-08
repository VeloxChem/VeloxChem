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


#include "SimdElectronRepulsionRecII.hpp"

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
#include "SimdTransformII.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ii_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ii_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(6535, nvalues);

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

            compute_prim_is_electron_repulsion_1(buffer, 305, 0, 140, 149, 221, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 317, 0, 149, 158, 233, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 329, 0, 158, 167, 245, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 341, 0, 167, 176, 257, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 353, 0, 176, 185, 269, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 365, 0, 185, 194, 281, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 377, 0, 194, 203, 293, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 389, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 392, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 395, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 398, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 401, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 404, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 407, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 410, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 413, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 416, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 419, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 422, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 425, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 428, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 431, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 434, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 437, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 440, 3, 17, 47, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 443, 3, 23, 59, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 446, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 449, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 452, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 461, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 470, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 479, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 488, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 491, 3, 47, 83, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 494, 3, 59, 92, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 497, 3, 62, 98, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 500, 3, 65, 104, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 509, 3, 68, 110, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 518, 3, 71, 116, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 527, 3, 74, 122, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 536, 3, 77, 128, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 545, 3, 80, 134, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 548, 3, 92, 158, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 551, 3, 98, 167, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 563, 3, 104, 176, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 575, 3, 110, 185, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 587, 3, 116, 194, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 599, 3, 122, 203, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 611, 3, 128, 212, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 620, 3, 158, 233, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 635, 3, 167, 245, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 650, 3, 176, 257, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 665, 3, 185, 269, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 680, 3, 194, 281, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 695, 3, 203, 293, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 710, 3, 233, 329, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 728, 3, 245, 341, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 746, 3, 257, 353, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 764, 3, 269, 365, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 782, 3, 281, 377, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 800, 3, 9, 10, 392, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 803, 3, 10, 11, 395, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 806, 3, 11, 12, 398, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 809, 3, 12, 13, 401, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 812, 3, 13, 14, 404, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 815, 3, 14, 15, 407, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 818, 3, 15, 16, 410, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 821, 3, 16, 17, 413, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 824, 0, 392, 803, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 827, 0, 395, 806, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 830, 0, 398, 809, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 833, 0, 401, 812, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 836, 0, 404, 815, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 839, 0, 407, 818, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 842, 0, 410, 821, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 845, 3, 416, 53, 56, 443, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 848, 3, 419, 56, 59, 446, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 851, 3, 422, 59, 62, 449, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 854, 3, 425, 62, 65, 452, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 857, 3, 428, 65, 68, 461, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 866, 3, 431, 68, 71, 470, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 875, 3, 434, 71, 74, 479, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 884, 3, 437, 74, 77, 488, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 887, 3, 440, 77, 80, 491, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 890, 0, 3, 446, 851, 86, 92, 497, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 899, 0, 3, 449, 854, 92, 98, 500, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 908, 0, 3, 452, 857, 98, 104, 509, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 923, 0, 3, 461, 866, 104, 110, 518, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 938, 0, 3, 470, 875, 110, 116, 527, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 953, 0, 3, 479, 884, 116, 122, 536, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 968, 0, 3, 488, 887, 122, 128, 545, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 977, 0, 3, 845, 848, 494, 890, 140, 149, 548, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 992, 0, 3, 848, 851, 497, 899, 149, 158, 551, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1007, 0, 3, 851, 854, 500, 908, 158, 167, 563, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 1031, 0, 3, 854, 857, 509, 923, 167, 176, 575, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1055, 0, 3, 857, 866, 518, 938, 176, 185, 587, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1079, 0, 3, 866, 875, 527, 953, 185, 194, 599, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_29(buffer, 1103, 0, 3, 875, 884, 536, 968, 194, 203, 611, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1124, 0, 3, 890, 899, 551, 1007, 221, 233, 635, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1157, 0, 3, 899, 908, 563, 1031, 233, 245, 650, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1190, 0, 3, 908, 923, 575, 1055, 245, 257, 665, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1223, 0, 3, 923, 938, 587, 1079, 257, 269, 680, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1256, 0, 3, 938, 953, 599, 1103, 269, 281, 695, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_2(buffer, 1289, 0, 3, 977, 992, 620, 1124, 305, 317, 710, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_5(buffer, 1325, 0, 3, 992, 1007, 635, 1157, 317, 329, 728, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1361, 0, 3, 1007, 1031, 650, 1190, 329, 341, 746, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1397, 0, 3, 1031, 1055, 665, 1223, 341, 353, 764, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_8(buffer, 1433, 0, 3, 1055, 1079, 680, 1256, 353, 365, 782, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1469, 3, 389, 392, 803, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1472, 3, 392, 395, 806, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1475, 3, 395, 398, 809, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1478, 3, 398, 401, 812, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1481, 3, 401, 404, 815, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1484, 3, 404, 407, 818, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1487, 3, 407, 410, 821, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1490, 0, 800, 1469, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1493, 0, 803, 1472, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1496, 0, 806, 1475, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1499, 0, 809, 1478, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1502, 0, 812, 1481, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1505, 0, 815, 1484, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1508, 0, 818, 1487, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1511, 3, 824, 443, 446, 851, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1514, 3, 827, 446, 449, 854, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_25(buffer, 1520, 3, 830, 449, 452, 857, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1526, 3, 833, 452, 461, 866, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1538, 3, 836, 461, 470, 875, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_27(buffer, 1550, 3, 839, 470, 479, 884, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_23(buffer, 1556, 3, 842, 479, 488, 887, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 1559, 0, 3, 851, 1514, 494, 497, 899, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_32(buffer, 1571, 0, 3, 854, 1520, 497, 500, 908, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_36(buffer, 1583, 0, 3, 857, 1526, 500, 509, 923, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_37(buffer, 1613, 0, 3, 866, 1538, 509, 518, 938, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_34(buffer, 1640, 0, 3, 875, 1550, 518, 527, 953, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_30(buffer, 1667, 0, 3, 884, 1556, 527, 536, 968, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_30(buffer, 1676, 0, 3, 1511, 1514, 899, 1571, 548, 551, 1007, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_34(buffer, 1694, 0, 3, 1514, 1520, 908, 1583, 551, 563, 1031, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_35(buffer, 1739, 0, 3, 1520, 1526, 923, 1613, 563, 575, 1055, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_36(buffer, 1790, 0, 3, 1526, 1538, 938, 1640, 575, 587, 1079, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_33(buffer, 1835, 0, 3, 1538, 1550, 953, 1667, 587, 599, 1103, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_16(buffer, 1868, 0, 3, 1559, 1571, 1007, 1694, 620, 635, 1157, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_17(buffer, 1928, 0, 3, 1571, 1583, 1031, 1739, 635, 650, 1190, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_18(buffer, 1991, 0, 3, 1583, 1613, 1055, 1790, 650, 665, 1223, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_19(buffer, 2057, 0, 3, 1613, 1640, 1079, 1835, 665, 680, 1256, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_5(buffer, 2117, 0, 3, 1676, 1694, 1157, 1928, 710, 728, 1361, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_6(buffer, 2189, 0, 3, 1694, 1739, 1190, 1991, 728, 746, 1397, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_7(buffer, 2261, 0, 3, 1739, 1790, 1223, 2057, 746, 764, 1433, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2333, 3, 800, 803, 1472, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2336, 3, 803, 806, 1475, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2339, 3, 806, 809, 1478, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2342, 3, 809, 812, 1481, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2345, 3, 812, 815, 1484, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2348, 3, 815, 818, 1487, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2351, 0, 1472, 2336, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2354, 0, 1475, 2339, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2357, 0, 1478, 2342, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2360, 0, 1481, 2345, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2363, 0, 1484, 2348, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2366, 3, 1490, 845, 848, 1511, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 2369, 3, 1493, 848, 851, 1514, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 2372, 3, 1496, 851, 854, 1520, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_32(buffer, 2381, 3, 1499, 854, 857, 1526, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_33(buffer, 2390, 3, 1502, 857, 866, 1538, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_34(buffer, 2405, 3, 1505, 866, 875, 1550, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_26(buffer, 2414, 3, 1508, 875, 884, 1556, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_34(buffer, 2417, 0, 3, 1514, 2372, 890, 899, 1571, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_35(buffer, 2426, 0, 3, 1520, 2381, 899, 908, 1583, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_36(buffer, 2441, 0, 3, 1526, 2390, 908, 923, 1613, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_37(buffer, 2486, 0, 3, 1538, 2405, 923, 938, 1640, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_33(buffer, 2525, 0, 3, 1550, 2414, 938, 953, 1667, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_29(buffer, 2534, 0, 3, 2366, 2369, 1559, 2417, 977, 992, 1676, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_30(buffer, 2555, 0, 3, 2369, 2372, 1571, 2426, 992, 1007, 1694, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_31(buffer, 2576, 0, 3, 2372, 2381, 1583, 2441, 1007, 1031, 1739, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_32(buffer, 2663, 0, 3, 2381, 2390, 1613, 2486, 1031, 1055, 1790, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_33(buffer, 2738, 0, 3, 2390, 2405, 1640, 2525, 1055, 1079, 1835, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_13(buffer, 2783, 0, 3, 2417, 2426, 1694, 2576, 1124, 1157, 1928, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_14(buffer, 2871, 0, 3, 2426, 2441, 1739, 2663, 1157, 1190, 1991, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_15(buffer, 2994, 0, 3, 2441, 2486, 1790, 2738, 1190, 1223, 2057, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_2(buffer, 3091, 0, 3, 2534, 2555, 1868, 2783, 1289, 1325, 2117, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_3(buffer, 3200, 0, 3, 2555, 2576, 1928, 2871, 1325, 1361, 2189, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_4(buffer, 3309, 0, 3, 2576, 2663, 1991, 2994, 1361, 1397, 2261, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3451, 3, 1469, 1472, 2336, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3454, 3, 1472, 1475, 2339, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3457, 3, 1475, 1478, 2342, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3460, 3, 1478, 1481, 2345, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3463, 3, 1481, 1484, 2348, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3466, 0, 2333, 3451, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3469, 0, 2336, 3454, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3472, 0, 2339, 3457, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3475, 0, 2342, 3460, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3478, 0, 2345, 3463, ncols, p);

            compute_prim_dh_electron_repulsion_32(buffer, 3481, 3, 2351, 1511, 1514, 2372, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_33(buffer, 3484, 3, 2354, 1514, 1520, 2381, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_34(buffer, 3493, 3, 2357, 1520, 1526, 2390, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_35(buffer, 3502, 3, 2360, 1526, 1538, 2405, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_31(buffer, 3511, 3, 2363, 1538, 1550, 2414, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_30(buffer, 3514, 0, 3, 2372, 3484, 1559, 1571, 2426, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_31(buffer, 3532, 0, 3, 2381, 3493, 1571, 1583, 2441, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_32(buffer, 3550, 0, 3, 2390, 3502, 1583, 1613, 2486, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_33(buffer, 3592, 0, 3, 2405, 3511, 1613, 1640, 2525, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_21(buffer, 3601, 0, 3, 3481, 3484, 2426, 3532, 1676, 1694, 2576, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_22(buffer, 3625, 0, 3, 3484, 3493, 2441, 3550, 1694, 1739, 2663, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_23(buffer, 3737, 0, 3, 3493, 3502, 2486, 3592, 1739, 1790, 2738, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_9(buffer, 3797, 0, 3, 3514, 3532, 2576, 3625, 1868, 1928, 2871, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_10(buffer, 4077, 0, 3, 3532, 3550, 2663, 3737, 1928, 1991, 2994, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_1(buffer, 4246, 0, 3, 3601, 3625, 2871, 4077, 2117, 2189, 3309, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_9(buffer, 4616, 3, 3466, 2366, 2369, 3481, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_25(buffer, 4619, 3, 3469, 2369, 2372, 3484, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_26(buffer, 4622, 3, 3472, 2372, 2381, 3493, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_27(buffer, 4625, 3, 3475, 2381, 2390, 3502, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_28(buffer, 4628, 3, 3478, 2390, 2405, 3511, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_20(buffer, 4631, 0, 3, 3484, 4622, 2417, 2426, 3532, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_21(buffer, 4640, 0, 3, 3493, 4625, 2426, 2441, 3550, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_22(buffer, 4649, 0, 3, 3502, 4628, 2441, 2486, 3592, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_11(buffer, 4658, 0, 3, 4616, 4619, 3514, 4631, 2534, 2555, 3601, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_12(buffer, 4682, 0, 3, 4619, 4622, 3532, 4640, 2555, 2576, 3625, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_13(buffer, 4706, 0, 3, 4622, 4625, 3550, 4649, 2576, 2663, 3737, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_5(buffer, 4766, 0, 3, 4631, 4640, 3625, 4706, 2783, 2871, 4077, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 4967, 0, 3, 4658, 4682, 3797, 4766, 3091, 3200, 4246, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5751, 4967, 784, ncols);
        }
    }

    simdtrf::transform_ii_tri(values, nvalues, buffer, 5751, nmax);
}

}  // namespace simdt2ceri
