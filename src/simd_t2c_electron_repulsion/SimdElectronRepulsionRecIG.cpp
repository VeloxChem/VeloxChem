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


#include "SimdElectronRepulsionRecIG.hpp"

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
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ig_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ig_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(9074, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 13, 14, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 14, 15, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 15, 16, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 99, 0, 18, 21, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 109, 0, 21, 24, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 24, 27, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 27, 30, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 30, 33, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 33, 36, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 36, 39, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 169, 0, 45, 51, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 184, 0, 51, 57, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 199, 0, 57, 63, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 214, 0, 63, 69, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 229, 0, 69, 75, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 244, 0, 75, 81, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 259, 0, 81, 87, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 274, 0, 99, 109, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 295, 0, 109, 119, 214, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 316, 0, 119, 129, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 337, 0, 129, 139, 244, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 358, 0, 139, 149, 259, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 379, 0, 169, 184, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 407, 0, 184, 199, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 435, 0, 199, 214, 316, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 463, 0, 214, 229, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 491, 0, 229, 244, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 519, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 522, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 525, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 528, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 531, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 534, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 537, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 540, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 543, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 552, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 561, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 570, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 579, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 588, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 597, 3, 15, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 606, 3, 16, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 615, 0, 3, 21, 552, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 633, 0, 3, 24, 561, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 651, 0, 3, 27, 570, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 669, 0, 3, 30, 579, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 687, 0, 3, 33, 588, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 705, 0, 3, 36, 597, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 723, 0, 3, 39, 606, 93, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 741, 0, 3, 57, 633, 109, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 771, 0, 3, 63, 651, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 801, 0, 3, 69, 669, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 831, 0, 3, 75, 687, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 861, 0, 3, 81, 705, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 891, 0, 3, 87, 723, 159, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 921, 0, 3, 109, 771, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 966, 0, 3, 119, 801, 214, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1011, 0, 3, 129, 831, 229, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1056, 0, 3, 139, 861, 244, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1101, 0, 3, 149, 891, 259, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1146, 0, 3, 199, 966, 295, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1209, 0, 3, 214, 1011, 316, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1272, 0, 3, 229, 1056, 337, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1335, 0, 3, 244, 1101, 358, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1398, 0, 3, 295, 1209, 435, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1482, 0, 3, 316, 1272, 463, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1566, 0, 3, 337, 1335, 491, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1650, 3, 9, 10, 522, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1656, 3, 10, 11, 525, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1662, 3, 11, 12, 528, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1668, 3, 12, 13, 531, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1674, 3, 13, 14, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1680, 3, 14, 15, 537, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1686, 3, 15, 16, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1692, 0, 3, 519, 1650, 552, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1710, 0, 3, 522, 1656, 561, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1728, 0, 3, 525, 1662, 570, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1746, 0, 3, 528, 1668, 579, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1764, 0, 3, 531, 1674, 588, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1782, 0, 3, 534, 1680, 597, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1800, 0, 3, 537, 1686, 606, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1818, 0, 3, 543, 1692, 45, 51, 615,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1854, 0, 3, 552, 1710, 51, 57, 633,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1890, 0, 3, 561, 1728, 57, 63, 651,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1926, 0, 3, 570, 1746, 63, 69, 669,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1962, 0, 3, 579, 1764, 69, 75, 687,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1998, 0, 3, 588, 1782, 75, 81, 705,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2034, 0, 3, 597, 1800, 81, 87, 723,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2070, 0, 3, 633, 1890, 99, 109, 771,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2130, 0, 3, 651, 1926, 109, 119, 801,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2190, 0, 3, 669, 1962, 119, 129, 831,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2250, 0, 3, 687, 1998, 129, 139, 861,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2310, 0, 3, 705, 2034, 139, 149, 891,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2370, 0, 3, 1818, 1854, 741, 2070, 169,
                                                 184, 921, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2460, 0, 3, 1854, 1890, 771, 2130, 184,
                                                 199, 966, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2550, 0, 3, 1890, 1926, 801, 2190, 199,
                                                 214, 1011, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2640, 0, 3, 1926, 1962, 831, 2250, 214,
                                                 229, 1056, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2730, 0, 3, 1962, 1998, 861, 2310, 229,
                                                 244, 1101, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2820, 0, 3, 2070, 2130, 966, 2550, 274,
                                                 295, 1209, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2946, 0, 3, 2130, 2190, 1011, 2640, 295,
                                                 316, 1272, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3072, 0, 3, 2190, 2250, 1056, 2730, 316,
                                                 337, 1335, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3198, 0, 3, 2370, 2460, 1146, 2820, 379,
                                                 407, 1398, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3366, 0, 3, 2460, 2550, 1209, 2946, 407,
                                                 435, 1482, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3534, 0, 3, 2550, 2640, 1272, 3072, 435,
                                                 463, 1566, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3702, 3, 519, 522, 1656, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3712, 3, 522, 525, 1662, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3722, 3, 525, 528, 1668, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3732, 3, 528, 531, 1674, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3742, 3, 531, 534, 1680, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3752, 3, 534, 537, 1686, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3762, 0, 3, 1650, 3702, 1710, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3792, 0, 3, 1656, 3712, 1728, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3822, 0, 3, 1662, 3722, 1746, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3852, 0, 3, 1668, 3732, 1764, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3882, 0, 3, 1674, 3742, 1782, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3912, 0, 3, 1680, 3752, 1800, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3942, 0, 3, 1710, 3792, 615, 633, 1890,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4002, 0, 3, 1728, 3822, 633, 651, 1926,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4062, 0, 3, 1746, 3852, 651, 669, 1962,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4122, 0, 3, 1764, 3882, 669, 687, 1998,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4182, 0, 3, 1782, 3912, 687, 705, 2034,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4242, 0, 3, 1890, 4002, 741, 771, 2130,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4342, 0, 3, 1926, 4062, 771, 801, 2190,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4442, 0, 3, 1962, 4122, 801, 831, 2250,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4542, 0, 3, 1998, 4182, 831, 861, 2310,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4642, 0, 3, 3942, 4002, 2130, 4342, 921,
                                                 966, 2550, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4792, 0, 3, 4002, 4062, 2190, 4442, 966,
                                                 1011, 2640, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4942, 0, 3, 4062, 4122, 2250, 4542,
                                                 1011, 1056, 2730, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5092, 0, 3, 4242, 4342, 2550, 4792,
                                                 1146, 1209, 2946, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5302, 0, 3, 4342, 4442, 2640, 4942,
                                                 1209, 1272, 3072, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 5512, 0, 3, 4642, 4792, 2946, 5302,
                                                 1398, 1482, 3534, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5792, 3, 1650, 1656, 3712, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5807, 3, 1656, 1662, 3722, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5822, 3, 1662, 1668, 3732, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5837, 3, 1668, 1674, 3742, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5852, 3, 1674, 1680, 3752, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 5867, 0, 3, 3702, 5792, 3792, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5912, 0, 3, 3712, 5807, 3822, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5957, 0, 3, 3722, 5822, 3852, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6002, 0, 3, 3732, 5837, 3882, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6047, 0, 3, 3742, 5852, 3912, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 6092, 0, 3, 3762, 5867, 1818, 1854,
                                                 3942, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6182, 0, 3, 3792, 5912, 1854, 1890,
                                                 4002, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6272, 0, 3, 3822, 5957, 1890, 1926,
                                                 4062, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6362, 0, 3, 3852, 6002, 1926, 1962,
                                                 4122, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6452, 0, 3, 3882, 6047, 1962, 1998,
                                                 4182, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6542, 0, 3, 4002, 6272, 2070, 2130,
                                                 4342, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6692, 0, 3, 4062, 6362, 2130, 2190,
                                                 4442, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6842, 0, 3, 4122, 6452, 2190, 2250,
                                                 4542, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 6992, 0, 3, 6092, 6182, 4242, 6542,
                                                 2370, 2460, 4642, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7217, 0, 3, 6182, 6272, 4342, 6692,
                                                 2460, 2550, 4792, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7442, 0, 3, 6272, 6362, 4442, 6842,
                                                 2550, 2640, 4942, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 7667, 0, 3, 6542, 6692, 4792, 7442,
                                                 2820, 2946, 5302, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 7982, 0, 3, 6992, 7217, 5092, 7667,
                                                 3198, 3366, 5512, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8402, 7982, 420, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 8822, 8402, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 8822, 9, nmax);
}

}  // namespace simdt2ceri
