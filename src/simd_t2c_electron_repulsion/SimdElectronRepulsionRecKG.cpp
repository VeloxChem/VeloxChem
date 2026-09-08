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


#include "SimdElectronRepulsionRecKG.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(13728, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 18, 21, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 21, 24, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 24, 27, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 27, 30, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 30, 33, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 33, 36, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 36, 39, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 39, 42, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 42, 45, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 195, 0, 51, 57, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 210, 0, 57, 63, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 225, 0, 63, 69, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 240, 0, 69, 75, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 255, 0, 75, 81, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 270, 0, 81, 87, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 285, 0, 87, 93, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 300, 0, 105, 115, 195, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 321, 0, 115, 125, 210, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 342, 0, 125, 135, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 363, 0, 135, 145, 240, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 384, 0, 145, 155, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 405, 0, 155, 165, 270, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 426, 0, 165, 175, 285, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 447, 0, 195, 210, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 475, 0, 210, 225, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 503, 0, 225, 240, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 531, 0, 240, 255, 405, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 559, 0, 255, 270, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 587, 0, 300, 321, 447, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 623, 0, 321, 342, 475, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 659, 0, 342, 363, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 695, 0, 363, 384, 531, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 731, 0, 384, 405, 559, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 767, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 770, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 773, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 776, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 779, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 782, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 785, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 788, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 791, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 800, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 809, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 818, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 827, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 836, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 845, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 854, 3, 16, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 863, 0, 3, 24, 791, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 881, 0, 3, 27, 800, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 899, 0, 3, 30, 809, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 917, 0, 3, 33, 818, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 935, 0, 3, 36, 827, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 953, 0, 3, 39, 836, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 971, 0, 3, 42, 845, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 989, 0, 3, 45, 854, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1007, 0, 3, 57, 881, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1037, 0, 3, 63, 899, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1067, 0, 3, 69, 917, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1097, 0, 3, 75, 935, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1127, 0, 3, 81, 953, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1157, 0, 3, 87, 971, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1187, 0, 3, 93, 989, 185, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1217, 0, 3, 125, 1037, 210, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1262, 0, 3, 135, 1067, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1307, 0, 3, 145, 1097, 240, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1352, 0, 3, 155, 1127, 255, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1397, 0, 3, 165, 1157, 270, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1442, 0, 3, 175, 1187, 285, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1487, 0, 3, 210, 1262, 342, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1550, 0, 3, 225, 1307, 363, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1613, 0, 3, 240, 1352, 384, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1676, 0, 3, 255, 1397, 405, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1739, 0, 3, 270, 1442, 426, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1802, 0, 3, 342, 1550, 475, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1886, 0, 3, 363, 1613, 503, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1970, 0, 3, 384, 1676, 531, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2054, 0, 3, 405, 1739, 559, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2138, 0, 3, 475, 1886, 659, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2246, 0, 3, 503, 1970, 695, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2354, 0, 3, 531, 2054, 731, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2462, 3, 9, 10, 770, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2468, 3, 10, 11, 773, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2474, 3, 11, 12, 776, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2480, 3, 12, 13, 779, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2486, 3, 13, 14, 782, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2492, 3, 14, 15, 785, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2498, 3, 15, 16, 788, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2504, 0, 3, 767, 2462, 800, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2522, 0, 3, 770, 2468, 809, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2540, 0, 3, 773, 2474, 818, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2558, 0, 3, 776, 2480, 827, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2576, 0, 3, 779, 2486, 836, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2594, 0, 3, 782, 2492, 845, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2612, 0, 3, 785, 2498, 854, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2630, 0, 3, 791, 2504, 51, 57, 881,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2666, 0, 3, 800, 2522, 57, 63, 899,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2702, 0, 3, 809, 2540, 63, 69, 917,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2738, 0, 3, 818, 2558, 69, 75, 935,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2774, 0, 3, 827, 2576, 75, 81, 953,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2810, 0, 3, 836, 2594, 81, 87, 971,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2846, 0, 3, 845, 2612, 87, 93, 989,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2882, 0, 3, 863, 2630, 105, 115, 1007,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2942, 0, 3, 881, 2666, 115, 125, 1037,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3002, 0, 3, 899, 2702, 125, 135, 1067,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3062, 0, 3, 917, 2738, 135, 145, 1097,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3122, 0, 3, 935, 2774, 145, 155, 1127,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3182, 0, 3, 953, 2810, 155, 165, 1157,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3242, 0, 3, 971, 2846, 165, 175, 1187,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3302, 0, 3, 2630, 2666, 1037, 3002, 195,
                                                 210, 1262, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3392, 0, 3, 2666, 2702, 1067, 3062, 210,
                                                 225, 1307, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3482, 0, 3, 2702, 2738, 1097, 3122, 225,
                                                 240, 1352, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3572, 0, 3, 2738, 2774, 1127, 3182, 240,
                                                 255, 1397, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3662, 0, 3, 2774, 2810, 1157, 3242, 255,
                                                 270, 1442, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3752, 0, 3, 2882, 2942, 1217, 3302, 300,
                                                 321, 1487, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3878, 0, 3, 2942, 3002, 1262, 3392, 321,
                                                 342, 1550, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4004, 0, 3, 3002, 3062, 1307, 3482, 342,
                                                 363, 1613, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4130, 0, 3, 3062, 3122, 1352, 3572, 363,
                                                 384, 1676, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4256, 0, 3, 3122, 3182, 1397, 3662, 384,
                                                 405, 1739, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4382, 0, 3, 3302, 3392, 1550, 4004, 447,
                                                 475, 1886, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4550, 0, 3, 3392, 3482, 1613, 4130, 475,
                                                 503, 1970, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4718, 0, 3, 3482, 3572, 1676, 4256, 503,
                                                 531, 2054, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 4886, 0, 3, 3752, 3878, 1802, 4382, 587,
                                                 623, 2138, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5102, 0, 3, 3878, 4004, 1886, 4550, 623,
                                                 659, 2246, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 5318, 0, 3, 4004, 4130, 1970, 4718, 659,
                                                 695, 2354, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5534, 3, 767, 770, 2468, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5544, 3, 770, 773, 2474, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5554, 3, 773, 776, 2480, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5564, 3, 776, 779, 2486, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5574, 3, 779, 782, 2492, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5584, 3, 782, 785, 2498, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 5594, 0, 3, 2462, 5534, 2522, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5624, 0, 3, 2468, 5544, 2540, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5654, 0, 3, 2474, 5554, 2558, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5684, 0, 3, 2480, 5564, 2576, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5714, 0, 3, 2486, 5574, 2594, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5744, 0, 3, 2492, 5584, 2612, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 5774, 0, 3, 2504, 5594, 863, 881, 2666,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5834, 0, 3, 2522, 5624, 881, 899, 2702,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5894, 0, 3, 2540, 5654, 899, 917, 2738,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5954, 0, 3, 2558, 5684, 917, 935, 2774,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6014, 0, 3, 2576, 5714, 935, 953, 2810,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6074, 0, 3, 2594, 5744, 953, 971, 2846,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6134, 0, 3, 2666, 5834, 1007, 1037,
                                                 3002, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6234, 0, 3, 2702, 5894, 1037, 1067,
                                                 3062, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6334, 0, 3, 2738, 5954, 1067, 1097,
                                                 3122, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6434, 0, 3, 2774, 6014, 1097, 1127,
                                                 3182, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6534, 0, 3, 2810, 6074, 1127, 1157,
                                                 3242, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6634, 0, 3, 5774, 5834, 3002, 6234,
                                                 1217, 1262, 3392, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6784, 0, 3, 5834, 5894, 3062, 6334,
                                                 1262, 1307, 3482, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6934, 0, 3, 5894, 5954, 3122, 6434,
                                                 1307, 1352, 3572, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7084, 0, 3, 5954, 6014, 3182, 6534,
                                                 1352, 1397, 3662, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7234, 0, 3, 6134, 6234, 3392, 6784,
                                                 1487, 1550, 4004, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7444, 0, 3, 6234, 6334, 3482, 6934,
                                                 1550, 1613, 4130, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7654, 0, 3, 6334, 6434, 3572, 7084,
                                                 1613, 1676, 4256, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 7864, 0, 3, 6634, 6784, 4004, 7444,
                                                 1802, 1886, 4550, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8144, 0, 3, 6784, 6934, 4130, 7654,
                                                 1886, 1970, 4718, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 8424, 0, 3, 7234, 7444, 4550, 8144,
                                                 2138, 2246, 5318, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8784, 3, 2462, 2468, 5544, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8799, 3, 2468, 2474, 5554, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8814, 3, 2474, 2480, 5564, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8829, 3, 2480, 2486, 5574, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8844, 3, 2486, 2492, 5584, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 8859, 0, 3, 5534, 8784, 5624, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8904, 0, 3, 5544, 8799, 5654, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8949, 0, 3, 5554, 8814, 5684, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8994, 0, 3, 5564, 8829, 5714, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9039, 0, 3, 5574, 8844, 5744, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 9084, 0, 3, 5594, 8859, 2630, 2666,
                                                 5834, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9174, 0, 3, 5624, 8904, 2666, 2702,
                                                 5894, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9264, 0, 3, 5654, 8949, 2702, 2738,
                                                 5954, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9354, 0, 3, 5684, 8994, 2738, 2774,
                                                 6014, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9444, 0, 3, 5714, 9039, 2774, 2810,
                                                 6074, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9534, 0, 3, 5774, 9084, 2882, 2942,
                                                 6134, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9684, 0, 3, 5834, 9174, 2942, 3002,
                                                 6234, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9834, 0, 3, 5894, 9264, 3002, 3062,
                                                 6334, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9984, 0, 3, 5954, 9354, 3062, 3122,
                                                 6434, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10134, 0, 3, 6014, 9444, 3122, 3182,
                                                 6534, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10284, 0, 3, 9084, 9174, 6234, 9834,
                                                 3302, 3392, 6784, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10509, 0, 3, 9174, 9264, 6334, 9984,
                                                 3392, 3482, 6934, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10734, 0, 3, 9264, 9354, 6434, 10134,
                                                 3482, 3572, 7084, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 10959, 0, 3, 9534, 9684, 6634, 10284,
                                                 3752, 3878, 7234, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 11274, 0, 3, 9684, 9834, 6784, 10509,
                                                 3878, 4004, 7444, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 11589, 0, 3, 9834, 9984, 6934, 10734,
                                                 4004, 4130, 7654, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 11904, 0, 3, 10284, 10509, 7444, 11589,
                                                 4382, 4550, 8144, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 12324, 0, 3, 10959, 11274, 7864, 11904,
                                                 4886, 5102, 8424, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 12864, 12324, 540, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 13404, 12864, 36, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 13404, 9, nmax);
}

}  // namespace simdt2ceri
