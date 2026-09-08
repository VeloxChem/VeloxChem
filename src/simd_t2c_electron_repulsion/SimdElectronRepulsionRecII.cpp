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
#include "SimdBoysFunc.hpp"

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
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ii_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 24691, 23543, 784, nvalues);

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

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 15, 16, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 16, 17, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 17, 18, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 20, 23, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 23, 26, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 26, 29, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 29, 32, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 32, 35, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 35, 38, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 38, 41, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 41, 44, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 44, 47, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 209, 0, 53, 59, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 224, 0, 59, 65, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 239, 0, 65, 71, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 254, 0, 71, 77, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 269, 0, 77, 83, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 284, 0, 83, 89, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 299, 0, 89, 95, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 314, 0, 95, 101, 189, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 329, 0, 101, 107, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 344, 0, 119, 129, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 365, 0, 129, 139, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 386, 0, 139, 149, 269, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 407, 0, 149, 159, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 428, 0, 159, 169, 299, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 449, 0, 169, 179, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 470, 0, 179, 189, 329, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 491, 0, 209, 224, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 519, 0, 224, 239, 365, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 547, 0, 239, 254, 386, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 575, 0, 254, 269, 407, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 603, 0, 269, 284, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 631, 0, 284, 299, 449, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 659, 0, 299, 314, 470, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 687, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 690, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 693, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 696, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 699, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 702, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 705, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 708, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 711, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 714, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 717, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 726, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 735, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 744, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 753, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 762, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 771, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 780, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 789, 3, 17, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 798, 3, 18, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 807, 0, 3, 23, 726, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 825, 0, 3, 26, 735, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 843, 0, 3, 29, 744, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 861, 0, 3, 32, 753, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 879, 0, 3, 35, 762, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 897, 0, 3, 38, 771, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 915, 0, 3, 41, 780, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 933, 0, 3, 44, 789, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 951, 0, 3, 47, 798, 113, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 969, 0, 3, 65, 825, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 999, 0, 3, 71, 843, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1029, 0, 3, 77, 861, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1059, 0, 3, 83, 879, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1089, 0, 3, 89, 897, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1119, 0, 3, 95, 915, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1149, 0, 3, 101, 933, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1179, 0, 3, 107, 951, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1209, 0, 3, 129, 999, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1254, 0, 3, 139, 1029, 254, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1299, 0, 3, 149, 1059, 269, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1344, 0, 3, 159, 1089, 284, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1389, 0, 3, 169, 1119, 299, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1434, 0, 3, 179, 1149, 314, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1479, 0, 3, 189, 1179, 329, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1524, 0, 3, 239, 1254, 365, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1587, 0, 3, 254, 1299, 386, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1650, 0, 3, 269, 1344, 407, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1713, 0, 3, 284, 1389, 428, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1776, 0, 3, 299, 1434, 449, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1839, 0, 3, 314, 1479, 470, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1902, 0, 3, 365, 1587, 547, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1986, 0, 3, 386, 1650, 575, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2070, 0, 3, 407, 1713, 603, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2154, 0, 3, 428, 1776, 631, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2238, 0, 3, 449, 1839, 659, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2322, 3, 9, 10, 690, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2328, 3, 10, 11, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2334, 3, 11, 12, 696, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2340, 3, 12, 13, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2346, 3, 13, 14, 702, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2352, 3, 14, 15, 705, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2358, 3, 15, 16, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2364, 3, 16, 17, 711, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2370, 3, 17, 18, 714, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2376, 0, 3, 687, 2322, 726, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2394, 0, 3, 690, 2328, 735, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2412, 0, 3, 693, 2334, 744, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2430, 0, 3, 696, 2340, 753, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2448, 0, 3, 699, 2346, 762, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2466, 0, 3, 702, 2352, 771, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2484, 0, 3, 705, 2358, 780, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2502, 0, 3, 708, 2364, 789, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2520, 0, 3, 711, 2370, 798, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2538, 0, 3, 717, 2376, 53, 59, 807,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2574, 0, 3, 726, 2394, 59, 65, 825,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2610, 0, 3, 735, 2412, 65, 71, 843,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2646, 0, 3, 744, 2430, 71, 77, 861,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2682, 0, 3, 753, 2448, 77, 83, 879,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2718, 0, 3, 762, 2466, 83, 89, 897,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2754, 0, 3, 771, 2484, 89, 95, 915,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2790, 0, 3, 780, 2502, 95, 101, 933,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2826, 0, 3, 789, 2520, 101, 107, 951,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2862, 0, 3, 825, 2610, 119, 129, 999,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2922, 0, 3, 843, 2646, 129, 139, 1029,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2982, 0, 3, 861, 2682, 139, 149, 1059,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3042, 0, 3, 879, 2718, 149, 159, 1089,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3102, 0, 3, 897, 2754, 159, 169, 1119,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3162, 0, 3, 915, 2790, 169, 179, 1149,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3222, 0, 3, 933, 2826, 179, 189, 1179,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3282, 0, 3, 2538, 2574, 969, 2862, 209,
                                                 224, 1209, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3372, 0, 3, 2574, 2610, 999, 2922, 224,
                                                 239, 1254, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3462, 0, 3, 2610, 2646, 1029, 2982, 239,
                                                 254, 1299, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3552, 0, 3, 2646, 2682, 1059, 3042, 254,
                                                 269, 1344, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3642, 0, 3, 2682, 2718, 1089, 3102, 269,
                                                 284, 1389, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3732, 0, 3, 2718, 2754, 1119, 3162, 284,
                                                 299, 1434, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3822, 0, 3, 2754, 2790, 1149, 3222, 299,
                                                 314, 1479, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3912, 0, 3, 2862, 2922, 1254, 3462, 344,
                                                 365, 1587, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4038, 0, 3, 2922, 2982, 1299, 3552, 365,
                                                 386, 1650, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4164, 0, 3, 2982, 3042, 1344, 3642, 386,
                                                 407, 1713, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4290, 0, 3, 3042, 3102, 1389, 3732, 407,
                                                 428, 1776, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4416, 0, 3, 3102, 3162, 1434, 3822, 428,
                                                 449, 1839, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4542, 0, 3, 3282, 3372, 1524, 3912, 491,
                                                 519, 1902, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4710, 0, 3, 3372, 3462, 1587, 4038, 519,
                                                 547, 1986, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4878, 0, 3, 3462, 3552, 1650, 4164, 547,
                                                 575, 2070, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5046, 0, 3, 3552, 3642, 1713, 4290, 575,
                                                 603, 2154, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5214, 0, 3, 3642, 3732, 1776, 4416, 603,
                                                 631, 2238, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5382, 3, 687, 690, 2328, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5392, 3, 690, 693, 2334, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5402, 3, 693, 696, 2340, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5412, 3, 696, 699, 2346, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5422, 3, 699, 702, 2352, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5432, 3, 702, 705, 2358, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5442, 3, 705, 708, 2364, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5452, 3, 708, 711, 2370, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 5462, 0, 3, 2322, 5382, 2394, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5492, 0, 3, 2328, 5392, 2412, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5522, 0, 3, 2334, 5402, 2430, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5552, 0, 3, 2340, 5412, 2448, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5582, 0, 3, 2346, 5422, 2466, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5612, 0, 3, 2352, 5432, 2484, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5642, 0, 3, 2358, 5442, 2502, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5672, 0, 3, 2364, 5452, 2520, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 5702, 0, 3, 2394, 5492, 807, 825, 2610,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5762, 0, 3, 2412, 5522, 825, 843, 2646,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5822, 0, 3, 2430, 5552, 843, 861, 2682,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5882, 0, 3, 2448, 5582, 861, 879, 2718,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5942, 0, 3, 2466, 5612, 879, 897, 2754,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6002, 0, 3, 2484, 5642, 897, 915, 2790,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6062, 0, 3, 2502, 5672, 915, 933, 2826,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6122, 0, 3, 2610, 5762, 969, 999, 2922,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6222, 0, 3, 2646, 5822, 999, 1029, 2982,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6322, 0, 3, 2682, 5882, 1029, 1059,
                                                 3042, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6422, 0, 3, 2718, 5942, 1059, 1089,
                                                 3102, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6522, 0, 3, 2754, 6002, 1089, 1119,
                                                 3162, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6622, 0, 3, 2790, 6062, 1119, 1149,
                                                 3222, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6722, 0, 3, 5702, 5762, 2922, 6222,
                                                 1209, 1254, 3462, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6872, 0, 3, 5762, 5822, 2982, 6322,
                                                 1254, 1299, 3552, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7022, 0, 3, 5822, 5882, 3042, 6422,
                                                 1299, 1344, 3642, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7172, 0, 3, 5882, 5942, 3102, 6522,
                                                 1344, 1389, 3732, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7322, 0, 3, 5942, 6002, 3162, 6622,
                                                 1389, 1434, 3822, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7472, 0, 3, 6122, 6222, 3462, 6872,
                                                 1524, 1587, 4038, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7682, 0, 3, 6222, 6322, 3552, 7022,
                                                 1587, 1650, 4164, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7892, 0, 3, 6322, 6422, 3642, 7172,
                                                 1650, 1713, 4290, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8102, 0, 3, 6422, 6522, 3732, 7322,
                                                 1713, 1776, 4416, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8312, 0, 3, 6722, 6872, 4038, 7682,
                                                 1902, 1986, 4878, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8592, 0, 3, 6872, 7022, 4164, 7892,
                                                 1986, 2070, 5046, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8872, 0, 3, 7022, 7172, 4290, 8102,
                                                 2070, 2154, 5214, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9152, 3, 2322, 2328, 5392, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9167, 3, 2328, 2334, 5402, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9182, 3, 2334, 2340, 5412, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9197, 3, 2340, 2346, 5422, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9212, 3, 2346, 2352, 5432, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9227, 3, 2352, 2358, 5442, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9242, 3, 2358, 2364, 5452, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 9257, 0, 3, 5382, 9152, 5492, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9302, 0, 3, 5392, 9167, 5522, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9347, 0, 3, 5402, 9182, 5552, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9392, 0, 3, 5412, 9197, 5582, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9437, 0, 3, 5422, 9212, 5612, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9482, 0, 3, 5432, 9227, 5642, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9527, 0, 3, 5442, 9242, 5672, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 9572, 0, 3, 5462, 9257, 2538, 2574,
                                                 5702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9662, 0, 3, 5492, 9302, 2574, 2610,
                                                 5762, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9752, 0, 3, 5522, 9347, 2610, 2646,
                                                 5822, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9842, 0, 3, 5552, 9392, 2646, 2682,
                                                 5882, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9932, 0, 3, 5582, 9437, 2682, 2718,
                                                 5942, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10022, 0, 3, 5612, 9482, 2718, 2754,
                                                 6002, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10112, 0, 3, 5642, 9527, 2754, 2790,
                                                 6062, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10202, 0, 3, 5762, 9752, 2862, 2922,
                                                 6222, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10352, 0, 3, 5822, 9842, 2922, 2982,
                                                 6322, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10502, 0, 3, 5882, 9932, 2982, 3042,
                                                 6422, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10652, 0, 3, 5942, 10022, 3042, 3102,
                                                 6522, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10802, 0, 3, 6002, 10112, 3102, 3162,
                                                 6622, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10952, 0, 3, 9572, 9662, 6122, 10202,
                                                 3282, 3372, 6722, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11177, 0, 3, 9662, 9752, 6222, 10352,
                                                 3372, 3462, 6872, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11402, 0, 3, 9752, 9842, 6322, 10502,
                                                 3462, 3552, 7022, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11627, 0, 3, 9842, 9932, 6422, 10652,
                                                 3552, 3642, 7172, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11852, 0, 3, 9932, 10022, 6522, 10802,
                                                 3642, 3732, 7322, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12077, 0, 3, 10202, 10352, 6872, 11402,
                                                 3912, 4038, 7682, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12392, 0, 3, 10352, 10502, 7022, 11627,
                                                 4038, 4164, 7892, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12707, 0, 3, 10502, 10652, 7172, 11852,
                                                 4164, 4290, 8102, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 13022, 0, 3, 10952, 11177, 7472, 12077,
                                                 4542, 4710, 8312, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 13442, 0, 3, 11177, 11402, 7682, 12392,
                                                 4710, 4878, 8592, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 13862, 0, 3, 11402, 11627, 7892, 12707,
                                                 4878, 5046, 8872, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14282, 3, 5382, 5392, 9167, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14303, 3, 5392, 5402, 9182, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14324, 3, 5402, 5412, 9197, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14345, 3, 5412, 5422, 9212, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14366, 3, 5422, 5432, 9227, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 14387, 3, 5432, 5442, 9242, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 14408, 0, 3, 9152, 14282, 9302, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 14471, 0, 3, 9167, 14303, 9347, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 14534, 0, 3, 9182, 14324, 9392, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 14597, 0, 3, 9197, 14345, 9437, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 14660, 0, 3, 9212, 14366, 9482, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 14723, 0, 3, 9227, 14387, 9527, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 14786, 0, 3, 9302, 14471, 5702, 5762,
                                                 9752, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14912, 0, 3, 9347, 14534, 5762, 5822,
                                                 9842, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 15038, 0, 3, 9392, 14597, 5822, 5882,
                                                 9932, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 15164, 0, 3, 9437, 14660, 5882, 5942,
                                                 10022, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 15290, 0, 3, 9482, 14723, 5942, 6002,
                                                 10112, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15416, 0, 3, 9752, 14912, 6122, 6222,
                                                 10352, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15626, 0, 3, 9842, 15038, 6222, 6322,
                                                 10502, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15836, 0, 3, 9932, 15164, 6322, 6422,
                                                 10652, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 16046, 0, 3, 10022, 15290, 6422, 6522,
                                                 10802, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16256, 0, 3, 14786, 14912, 10352, 15626,
                                                 6722, 6872, 11402, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16571, 0, 3, 14912, 15038, 10502, 15836,
                                                 6872, 7022, 11627, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16886, 0, 3, 15038, 15164, 10652, 16046,
                                                 7022, 7172, 11852, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 17201, 0, 3, 15416, 15626, 11402, 16571,
                                                 7472, 7682, 12392, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 17642, 0, 3, 15626, 15836, 11627, 16886,
                                                 7682, 7892, 12707, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 18083, 0, 3, 16256, 16571, 12392, 17642,
                                                 8312, 8592, 13862, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18671, 3, 9152, 9167, 14303, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18699, 3, 9167, 9182, 14324, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18727, 3, 9182, 9197, 14345, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18755, 3, 9197, 9212, 14366, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18783, 3, 9212, 9227, 14387, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 18811, 0, 3, 14282, 18671, 14471, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18895, 0, 3, 14303, 18699, 14534, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18979, 0, 3, 14324, 18727, 14597, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 19063, 0, 3, 14345, 18755, 14660, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 19147, 0, 3, 14366, 18783, 14723, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 19231, 0, 3, 14408, 18811, 9572, 9662,
                                                 14786, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19399, 0, 3, 14471, 18895, 9662, 9752,
                                                 14912, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19567, 0, 3, 14534, 18979, 9752, 9842,
                                                 15038, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19735, 0, 3, 14597, 19063, 9842, 9932,
                                                 15164, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19903, 0, 3, 14660, 19147, 9932, 10022,
                                                 15290, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20071, 0, 3, 14912, 19567, 10202, 10352,
                                                 15626, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20351, 0, 3, 15038, 19735, 10352, 10502,
                                                 15836, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20631, 0, 3, 15164, 19903, 10502, 10652,
                                                 16046, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 20911, 0, 3, 19231, 19399, 15416, 20071,
                                                 10952, 11177, 16256, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 21331, 0, 3, 19399, 19567, 15626, 20351,
                                                 11177, 11402, 16571, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 21751, 0, 3, 19567, 19735, 15836, 20631,
                                                 11402, 11627, 16886, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 22171, 0, 3, 20071, 20351, 16571, 21751,
                                                 12077, 12392, 17642, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 22759, 0, 3, 20911, 21331, 17201, 22171,
                                                 13022, 13442, 18083, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 23543, 22759, 784, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 24327, 23543, 28, nmax);

    simdtrf::transform_i_outer_tri(values, nvalues, buffer, 24327, nmax);
}

}  // namespace simdt2ceri
