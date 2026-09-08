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


#include "SimdElectronRepulsionRecKF.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformKF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kf_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2262, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 47, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 50, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 53, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 12, 13, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 13, 14, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 14, 15, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 68, 0, 17, 20, 47, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 74, 0, 20, 23, 50, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 80, 0, 23, 26, 53, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_3(buffer, 86, 0, 26, 29, 56, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 94, 0, 29, 32, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 32, 35, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 106, 0, 35, 38, 65, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 112, 0, 44, 47, 74, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 121, 0, 47, 50, 80, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_10(buffer, 130, 0, 50, 53, 86, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_4(buffer, 141, 0, 53, 56, 94, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 152, 0, 56, 59, 100, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 161, 0, 59, 62, 106, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 170, 0, 68, 74, 121, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_16(buffer, 182, 0, 74, 80, 130, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_17(buffer, 194, 0, 80, 86, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_13(buffer, 211, 0, 86, 94, 152, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 226, 0, 94, 100, 161, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 239, 0, 112, 121, 182, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_12(buffer, 251, 0, 121, 130, 194, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_13(buffer, 266, 0, 130, 141, 211, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_14(buffer, 287, 0, 141, 152, 226, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 307, 0, 170, 182, 251, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_4(buffer, 322, 0, 182, 194, 266, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_5(buffer, 337, 0, 194, 211, 287, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 361, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 364, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 367, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 370, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 373, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 376, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 379, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 382, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 385, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 388, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 391, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 394, 3, 14, 38, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 397, 3, 20, 47, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 400, 3, 23, 50, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 403, 3, 26, 53, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 406, 3, 29, 56, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 409, 3, 32, 59, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 412, 3, 35, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 415, 3, 38, 65, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 418, 3, 44, 68, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 421, 3, 47, 74, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 424, 3, 50, 80, ncols, p);

            compute_prim_fp_electron_repulsion_14(buffer, 427, 3, 53, 86, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 430, 3, 56, 94, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 433, 3, 59, 100, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 436, 3, 62, 106, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 439, 3, 74, 121, ncols, p);

            compute_prim_gp_electron_repulsion_16(buffer, 442, 3, 80, 130, ncols, p);

            compute_prim_gp_electron_repulsion_17(buffer, 445, 0, 3, 86, 430, 141, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 459, 3, 94, 152, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 468, 3, 100, 161, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 471, 3, 112, 170, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 474, 3, 121, 182, ncols, p);

            compute_prim_hp_electron_repulsion_16(buffer, 477, 0, 3, 130, 445, 194, ncols, p);

            compute_prim_hp_electron_repulsion_10(buffer, 502, 0, 3, 141, 459, 211, ncols, p);

            compute_prim_hp_electron_repulsion_14(buffer, 524, 3, 152, 226, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 533, 3, 182, 251, ncols, p);

            compute_prim_ip_electron_repulsion_9(buffer, 551, 0, 3, 194, 502, 266, ncols, p);

            compute_prim_ip_electron_repulsion_10(buffer, 590, 0, 3, 211, 524, 287, ncols, p);

            compute_prim_kp_electron_repulsion_2(buffer, 617, 3, 239, 307, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 638, 3, 251, 322, ncols, p);

            compute_prim_kp_electron_repulsion_4(buffer, 659, 0, 3, 266, 590, 337, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 708, 3, 9, 10, 364, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 711, 3, 10, 11, 367, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 714, 3, 11, 12, 370, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 717, 3, 12, 13, 373, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 720, 3, 13, 14, 376, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 723, 0, 361, 708, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 726, 0, 364, 711, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 729, 0, 367, 714, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 732, 0, 370, 717, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 735, 0, 373, 720, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 738, 3, 379, 44, 47, 400, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 741, 3, 382, 47, 50, 403, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 744, 3, 385, 50, 53, 406, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 747, 3, 388, 53, 56, 409, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 750, 3, 391, 56, 59, 412, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 753, 3, 394, 59, 62, 415, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 756, 0, 3, 400, 741, 68, 74, 424, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 765, 0, 3, 403, 744, 74, 80, 427, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_25(buffer, 774, 0, 3, 406, 747, 80, 86, 430, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_24(buffer, 783, 0, 3, 409, 750, 86, 94, 433, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 792, 0, 3, 412, 753, 94, 100, 436, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 801, 0, 3, 738, 741, 424, 765, 112, 121, 442, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_33(buffer, 816, 0, 3, 741, 744, 427, 774, 121, 130, 445, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_34(buffer, 831, 0, 3, 744, 747, 430, 783, 130, 141, 459, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_35(buffer, 852, 0, 3, 747, 750, 433, 792, 141, 152, 468, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_25(buffer, 867, 0, 3, 756, 765, 442, 816, 170, 182, 477, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_26(buffer, 891, 0, 3, 765, 774, 445, 831, 182, 194, 502, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_27(buffer, 937, 0, 3, 774, 783, 459, 852, 194, 211, 524, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_10(buffer, 967, 0, 3, 801, 816, 477, 891, 239, 251, 551, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_11(buffer, 1063, 0, 3, 816, 831, 502, 937, 251, 266, 590, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_1(buffer, 1125, 0, 3, 867, 891, 551, 1063, 307, 322, 659, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1266, 3, 723, 397, 400, 741, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1269, 3, 726, 400, 403, 744, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1272, 3, 729, 403, 406, 747, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1275, 3, 732, 406, 409, 750, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 1278, 3, 735, 409, 412, 753, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1281, 0, 3, 738, 1266, 418, 421, 756, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1290, 0, 3, 741, 1269, 421, 424, 765, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1299, 0, 3, 744, 1272, 424, 427, 774, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1308, 0, 3, 747, 1275, 427, 430, 783, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 1317, 0, 3, 750, 1278, 430, 433, 792, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 1326, 0, 3, 1266, 1269, 765, 1299, 439, 442, 816, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_39(buffer, 1341, 0, 3, 1269, 1272, 774, 1308, 442, 445, 831, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_40(buffer, 1356, 0, 3, 1272, 1275, 783, 1317, 445, 459, 852, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_23(buffer, 1371, 0, 3, 1281, 1290, 801, 1326, 471, 474, 867, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_24(buffer, 1395, 0, 3, 1290, 1299, 816, 1341, 474, 477, 891, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_25(buffer, 1419, 0, 3, 1299, 1308, 831, 1356, 477, 502, 937, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_10(buffer, 1449, 0, 3, 1326, 1341, 891, 1419, 533, 551, 1063, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 1542, 0, 3, 1371, 1395, 967, 1449, 617, 638, 1125, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1902, 1542, 360, ncols);
        }
    }

    simdtrf::transform_kf(values, nvalues, buffer, 1902, nmax);
}

}  // namespace simdt2ceri
