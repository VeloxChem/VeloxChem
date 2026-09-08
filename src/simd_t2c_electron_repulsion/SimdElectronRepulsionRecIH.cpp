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


#include "SimdElectronRepulsionRecIH.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformIH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ih_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ih_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(4310, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 7, 8, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 8, 9, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 9, 10, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 10, 11, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 11, 12, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 12, 13, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 13, 14, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 14, 15, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 15, 16, 45, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 75, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 81, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 87, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 93, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 99, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 105, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 111, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 117, 0, 39, 42, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 123, 0, 48, 51, 81, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 132, 0, 51, 54, 87, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 141, 0, 54, 57, 93, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 150, 0, 57, 60, 99, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 159, 0, 60, 63, 105, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 168, 0, 63, 66, 111, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 177, 0, 66, 69, 117, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 186, 0, 75, 81, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 195, 0, 81, 87, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 207, 0, 87, 93, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 219, 0, 93, 99, 159, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 231, 0, 99, 105, 168, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 243, 0, 105, 111, 177, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 255, 0, 123, 132, 195, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 267, 0, 132, 141, 207, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 279, 0, 141, 150, 219, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 291, 0, 150, 159, 231, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 303, 0, 159, 168, 243, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 315, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 318, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 321, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 324, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 327, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 330, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 333, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 336, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 339, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 342, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 345, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 348, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 351, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 354, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 357, 3, 15, 42, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 360, 3, 18, 48, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 363, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 366, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 369, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 372, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 381, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 390, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 399, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 402, 3, 42, 72, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 405, 3, 51, 81, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 408, 3, 54, 87, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 411, 3, 57, 93, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 420, 3, 60, 99, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 429, 3, 63, 105, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 438, 3, 66, 111, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 447, 3, 69, 117, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 450, 3, 75, 123, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 453, 3, 81, 132, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 456, 3, 87, 141, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 468, 3, 93, 150, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 480, 3, 99, 159, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 492, 3, 105, 168, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 504, 3, 111, 177, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 513, 3, 132, 195, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 528, 3, 141, 207, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 543, 3, 150, 219, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 558, 3, 159, 231, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 573, 3, 168, 243, ncols, p);

            compute_prim_ip_electron_repulsion_2(buffer, 588, 3, 186, 255, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 606, 3, 195, 267, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 624, 3, 207, 279, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 642, 3, 219, 291, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 660, 3, 231, 303, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 678, 3, 8, 9, 318, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 681, 3, 9, 10, 321, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 684, 3, 10, 11, 324, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 687, 3, 11, 12, 327, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 690, 3, 12, 13, 330, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 693, 3, 13, 14, 333, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 696, 3, 14, 15, 336, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 699, 0, 315, 678, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 702, 0, 318, 681, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 705, 0, 321, 684, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 708, 0, 324, 687, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 711, 0, 327, 690, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 714, 0, 330, 693, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 717, 0, 333, 696, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 720, 3, 339, 48, 51, 366, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 723, 3, 342, 51, 54, 369, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 726, 3, 345, 54, 57, 372, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 729, 3, 348, 57, 60, 381, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 738, 3, 351, 60, 63, 390, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 747, 3, 354, 63, 66, 399, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 750, 3, 357, 66, 69, 402, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 753, 0, 3, 366, 723, 75, 81, 408, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 762, 0, 3, 369, 726, 81, 87, 411, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 771, 0, 3, 372, 729, 87, 93, 420, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 786, 0, 3, 381, 738, 93, 99, 429, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 801, 0, 3, 390, 747, 99, 105, 438, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 816, 0, 3, 399, 750, 105, 111, 447, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 825, 0, 3, 720, 723, 408, 762, 123, 132, 456, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 840, 0, 3, 723, 726, 411, 771, 132, 141, 468, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_30(buffer, 864, 0, 3, 726, 729, 420, 786, 141, 150, 480, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 891, 0, 3, 729, 738, 429, 801, 150, 159, 492, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_29(buffer, 915, 0, 3, 738, 747, 438, 816, 159, 168, 504, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_15(buffer, 936, 0, 3, 753, 762, 456, 840, 186, 195, 528, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_16(buffer, 966, 0, 3, 762, 771, 468, 864, 195, 207, 543, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_17(buffer, 999, 0, 3, 771, 786, 480, 891, 207, 219, 558, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1035, 0, 3, 786, 801, 492, 915, 219, 231, 573, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_5(buffer, 1068, 0, 3, 825, 840, 528, 966, 255, 267, 624, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_6(buffer, 1104, 0, 3, 840, 864, 543, 999, 267, 279, 642, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_7(buffer, 1140, 0, 3, 864, 891, 558, 1035, 279, 291, 660, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1176, 3, 315, 318, 681, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1179, 3, 318, 321, 684, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1182, 3, 321, 324, 687, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1185, 3, 324, 327, 690, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1188, 3, 327, 330, 693, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1191, 3, 330, 333, 696, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1194, 0, 681, 1179, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1197, 0, 684, 1182, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1200, 0, 687, 1185, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1203, 0, 690, 1188, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1206, 0, 693, 1191, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1209, 3, 699, 360, 363, 720, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1212, 3, 702, 363, 366, 723, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 1215, 3, 705, 366, 369, 726, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_25(buffer, 1221, 3, 708, 369, 372, 729, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_26(buffer, 1227, 3, 711, 372, 381, 738, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_27(buffer, 1239, 3, 714, 381, 390, 747, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_23(buffer, 1245, 3, 717, 390, 399, 750, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_31(buffer, 1248, 0, 3, 723, 1215, 405, 408, 762, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_32(buffer, 1257, 0, 3, 726, 1221, 408, 411, 771, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_33(buffer, 1269, 0, 3, 729, 1227, 411, 420, 786, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_34(buffer, 1302, 0, 3, 738, 1239, 420, 429, 801, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_30(buffer, 1329, 0, 3, 747, 1245, 429, 438, 816, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_29(buffer, 1338, 0, 3, 1209, 1212, 753, 1248, 450, 453, 825, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_30(buffer, 1356, 0, 3, 1212, 1215, 762, 1257, 453, 456, 840, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_31(buffer, 1374, 0, 3, 1215, 1221, 771, 1269, 456, 468, 864, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_32(buffer, 1437, 0, 3, 1221, 1227, 786, 1302, 468, 480, 891, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_33(buffer, 1491, 0, 3, 1227, 1239, 801, 1329, 480, 492, 915, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_13(buffer, 1524, 0, 3, 1248, 1257, 840, 1374, 513, 528, 966, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_14(buffer, 1584, 0, 3, 1257, 1269, 864, 1437, 528, 543, 999, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_15(buffer, 1671, 0, 3, 1269, 1302, 891, 1491, 543, 558, 1035, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_2(buffer, 1740, 0, 3, 1338, 1356, 936, 1524, 588, 606, 1068, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_3(buffer, 1812, 0, 3, 1356, 1374, 966, 1584, 606, 624, 1104, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_4(buffer, 1884, 0, 3, 1374, 1437, 999, 1671, 624, 642, 1140, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1989, 3, 678, 681, 1179, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1992, 3, 681, 684, 1182, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1995, 3, 684, 687, 1185, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1998, 3, 687, 690, 1188, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2001, 3, 690, 693, 1191, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2004, 0, 1176, 1989, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2007, 0, 1179, 1992, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2010, 0, 1182, 1995, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2013, 0, 1185, 1998, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2016, 0, 1188, 2001, ncols, p);

            compute_prim_dg_electron_repulsion_27(buffer, 2019, 3, 1194, 720, 723, 1215, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_28(buffer, 2022, 3, 1197, 723, 726, 1221, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_29(buffer, 2028, 3, 1200, 726, 729, 1227, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_30(buffer, 2034, 3, 1203, 729, 738, 1239, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_26(buffer, 2040, 3, 1206, 738, 747, 1245, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_30(buffer, 2043, 0, 3, 1215, 2022, 753, 762, 1257, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_31(buffer, 2058, 0, 3, 1221, 2028, 762, 771, 1269, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_32(buffer, 2073, 0, 3, 1227, 2034, 771, 786, 1302, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_33(buffer, 2103, 0, 3, 1239, 2040, 786, 801, 1329, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_26(buffer, 2112, 0, 3, 2019, 2022, 1257, 2058, 825, 840, 1374, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_27(buffer, 2133, 0, 3, 2022, 2028, 1269, 2073, 840, 864, 1437, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_28(buffer, 2215, 0, 3, 2028, 2034, 1302, 2103, 864, 891, 1491, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_11(buffer, 2260, 0, 3, 2043, 2058, 1374, 2133, 936, 966, 1584, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_12(buffer, 2461, 0, 3, 2058, 2073, 1437, 2215, 966, 999, 1671, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_1(buffer, 2582, 0, 3, 2112, 2133, 1584, 2461, 1068, 1104, 1884, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_9(buffer, 2855, 3, 2004, 1209, 1212, 2019, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_28(buffer, 2858, 3, 2007, 1212, 1215, 2022, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_29(buffer, 2861, 3, 2010, 1215, 1221, 2028, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_30(buffer, 2864, 3, 2013, 1221, 1227, 2034, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_31(buffer, 2867, 3, 2016, 1227, 1239, 2040, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_27(buffer, 2870, 0, 3, 2022, 2861, 1248, 1257, 2058, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_28(buffer, 2879, 0, 3, 2028, 2864, 1257, 1269, 2073, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_29(buffer, 2888, 0, 3, 2034, 2867, 1269, 1302, 2103, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_18(buffer, 2897, 0, 3, 2855, 2858, 2043, 2870, 1338, 1356, 2112, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_19(buffer, 2918, 0, 3, 2858, 2861, 2058, 2879, 1356, 1374, 2133, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_20(buffer, 2939, 0, 3, 2861, 2864, 2073, 2888, 1374, 1437, 2215, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_8(buffer, 2984, 0, 3, 2870, 2879, 2133, 2939, 1524, 1584, 2461, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 3134, 0, 3, 2897, 2918, 2260, 2984, 1740, 1812, 2582, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3722, 3134, 588, ncols);
        }
    }

    simdtrf::transform_ih(values, nvalues, buffer, 3722, nmax);
}

}  // namespace simdt2ceri
