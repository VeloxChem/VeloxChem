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
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

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

    auto buffer = CSimdMatrix(15450, nvalues);

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

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 7, 8, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 8, 9, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 9, 10, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 10, 11, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 11, 12, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 12, 13, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 13, 14, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 14, 15, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 15, 16, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 18, 21, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 21, 24, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 24, 27, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 27, 30, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 30, 33, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 33, 36, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 36, 39, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 39, 42, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 182, 0, 48, 54, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 197, 0, 54, 60, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 212, 0, 60, 66, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 227, 0, 66, 72, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 242, 0, 72, 78, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 257, 0, 78, 84, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 272, 0, 84, 90, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 287, 0, 102, 112, 197, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 308, 0, 112, 122, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 329, 0, 122, 132, 227, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 350, 0, 132, 142, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 371, 0, 142, 152, 257, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 392, 0, 152, 162, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 413, 0, 182, 197, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 441, 0, 197, 212, 329, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 469, 0, 212, 227, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 497, 0, 227, 242, 371, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 242, 257, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 553, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 556, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 559, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 562, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 565, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 568, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 571, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 574, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 577, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 580, 3, 8, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 589, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 598, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 607, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 616, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 625, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 634, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 643, 3, 15, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 652, 3, 16, 45, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 661, 0, 3, 18, 580, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 679, 0, 3, 21, 589, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 697, 0, 3, 24, 598, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 715, 0, 3, 27, 607, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 733, 0, 3, 30, 616, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 751, 0, 3, 33, 625, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 769, 0, 3, 36, 634, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 787, 0, 3, 39, 643, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 805, 0, 3, 42, 652, 96, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 823, 0, 3, 54, 697, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 853, 0, 3, 60, 715, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 883, 0, 3, 66, 733, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 913, 0, 3, 72, 751, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 943, 0, 3, 78, 769, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 973, 0, 3, 84, 787, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1003, 0, 3, 90, 805, 172, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1033, 0, 3, 102, 823, 182, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1078, 0, 3, 112, 853, 197, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1123, 0, 3, 122, 883, 212, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1168, 0, 3, 132, 913, 227, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1213, 0, 3, 142, 943, 242, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1258, 0, 3, 152, 973, 257, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1303, 0, 3, 162, 1003, 272, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1348, 0, 3, 197, 1123, 308, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1411, 0, 3, 212, 1168, 329, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1474, 0, 3, 227, 1213, 350, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1537, 0, 3, 242, 1258, 371, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1600, 0, 3, 257, 1303, 392, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1663, 0, 3, 287, 1348, 413, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1747, 0, 3, 308, 1411, 441, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1831, 0, 3, 329, 1474, 469, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1915, 0, 3, 350, 1537, 497, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1999, 0, 3, 371, 1600, 525, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2083, 3, 8, 9, 556, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2089, 3, 9, 10, 559, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2095, 3, 10, 11, 562, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2101, 3, 11, 12, 565, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2107, 3, 12, 13, 568, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2113, 3, 13, 14, 571, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2119, 3, 14, 15, 574, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2125, 3, 15, 16, 577, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2131, 0, 3, 553, 2083, 589, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2149, 0, 3, 556, 2089, 598, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2167, 0, 3, 559, 2095, 607, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2185, 0, 3, 562, 2101, 616, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2203, 0, 3, 565, 2107, 625, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2221, 0, 3, 568, 2113, 634, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2239, 0, 3, 571, 2119, 643, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2257, 0, 3, 574, 2125, 652, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2275, 0, 3, 589, 2149, 48, 54, 697,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2311, 0, 3, 598, 2167, 54, 60, 715,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2347, 0, 3, 607, 2185, 60, 66, 733,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2383, 0, 3, 616, 2203, 66, 72, 751,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2419, 0, 3, 625, 2221, 72, 78, 769,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2455, 0, 3, 634, 2239, 78, 84, 787,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2491, 0, 3, 643, 2257, 84, 90, 805,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2527, 0, 3, 697, 2311, 102, 112, 853,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2587, 0, 3, 715, 2347, 112, 122, 883,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2647, 0, 3, 733, 2383, 122, 132, 913,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2707, 0, 3, 751, 2419, 132, 142, 943,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2767, 0, 3, 769, 2455, 142, 152, 973,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2827, 0, 3, 787, 2491, 152, 162, 1003,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2887, 0, 3, 2275, 2311, 853, 2587, 182,
                                                 197, 1123, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2977, 0, 3, 2311, 2347, 883, 2647, 197,
                                                 212, 1168, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3067, 0, 3, 2347, 2383, 913, 2707, 212,
                                                 227, 1213, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3157, 0, 3, 2383, 2419, 943, 2767, 227,
                                                 242, 1258, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3247, 0, 3, 2419, 2455, 973, 2827, 242,
                                                 257, 1303, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3337, 0, 3, 2527, 2587, 1123, 2977, 287,
                                                 308, 1411, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3463, 0, 3, 2587, 2647, 1168, 3067, 308,
                                                 329, 1474, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3589, 0, 3, 2647, 2707, 1213, 3157, 329,
                                                 350, 1537, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3715, 0, 3, 2707, 2767, 1258, 3247, 350,
                                                 371, 1600, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3841, 0, 3, 2887, 2977, 1411, 3463, 413,
                                                 441, 1831, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4009, 0, 3, 2977, 3067, 1474, 3589, 441,
                                                 469, 1915, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4177, 0, 3, 3067, 3157, 1537, 3715, 469,
                                                 497, 1999, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4345, 3, 553, 556, 2089, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4355, 3, 556, 559, 2095, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4365, 3, 559, 562, 2101, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4375, 3, 562, 565, 2107, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4385, 3, 565, 568, 2113, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4395, 3, 568, 571, 2119, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4405, 3, 571, 574, 2125, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4415, 0, 3, 2083, 4345, 2149, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4445, 0, 3, 2089, 4355, 2167, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4475, 0, 3, 2095, 4365, 2185, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4505, 0, 3, 2101, 4375, 2203, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4535, 0, 3, 2107, 4385, 2221, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4565, 0, 3, 2113, 4395, 2239, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4595, 0, 3, 2119, 4405, 2257, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4625, 0, 3, 2131, 4415, 661, 679, 2275,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4685, 0, 3, 2149, 4445, 679, 697, 2311,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4745, 0, 3, 2167, 4475, 697, 715, 2347,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4805, 0, 3, 2185, 4505, 715, 733, 2383,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4865, 0, 3, 2203, 4535, 733, 751, 2419,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4925, 0, 3, 2221, 4565, 751, 769, 2455,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4985, 0, 3, 2239, 4595, 769, 787, 2491,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5045, 0, 3, 2311, 4745, 823, 853, 2587,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5145, 0, 3, 2347, 4805, 853, 883, 2647,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5245, 0, 3, 2383, 4865, 883, 913, 2707,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5345, 0, 3, 2419, 4925, 913, 943, 2767,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5445, 0, 3, 2455, 4985, 943, 973, 2827,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5545, 0, 3, 4625, 4685, 2527, 5045,
                                                 1033, 1078, 2887, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5695, 0, 3, 4685, 4745, 2587, 5145,
                                                 1078, 1123, 2977, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5845, 0, 3, 4745, 4805, 2647, 5245,
                                                 1123, 1168, 3067, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5995, 0, 3, 4805, 4865, 2707, 5345,
                                                 1168, 1213, 3157, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6145, 0, 3, 4865, 4925, 2767, 5445,
                                                 1213, 1258, 3247, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6295, 0, 3, 5045, 5145, 2977, 5845,
                                                 1348, 1411, 3463, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6505, 0, 3, 5145, 5245, 3067, 5995,
                                                 1411, 1474, 3589, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6715, 0, 3, 5245, 5345, 3157, 6145,
                                                 1474, 1537, 3715, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 6925, 0, 3, 5545, 5695, 3337, 6295,
                                                 1663, 1747, 3841, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 7205, 0, 3, 5695, 5845, 3463, 6505,
                                                 1747, 1831, 4009, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 7485, 0, 3, 5845, 5995, 3589, 6715,
                                                 1831, 1915, 4177, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7765, 3, 2083, 2089, 4355, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7780, 3, 2089, 2095, 4365, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7795, 3, 2095, 2101, 4375, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7810, 3, 2101, 2107, 4385, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7825, 3, 2107, 2113, 4395, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7840, 3, 2113, 2119, 4405, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 7855, 0, 3, 4345, 7765, 4445, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7900, 0, 3, 4355, 7780, 4475, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7945, 0, 3, 4365, 7795, 4505, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7990, 0, 3, 4375, 7810, 4535, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8035, 0, 3, 4385, 7825, 4565, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8080, 0, 3, 4395, 7840, 4595, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 8125, 0, 3, 4445, 7900, 2275, 2311,
                                                 4745, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8215, 0, 3, 4475, 7945, 2311, 2347,
                                                 4805, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8305, 0, 3, 4505, 7990, 2347, 2383,
                                                 4865, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8395, 0, 3, 4535, 8035, 2383, 2419,
                                                 4925, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8485, 0, 3, 4565, 8080, 2419, 2455,
                                                 4985, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8575, 0, 3, 4745, 8215, 2527, 2587,
                                                 5145, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8725, 0, 3, 4805, 8305, 2587, 2647,
                                                 5245, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8875, 0, 3, 4865, 8395, 2647, 2707,
                                                 5345, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9025, 0, 3, 4925, 8485, 2707, 2767,
                                                 5445, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9175, 0, 3, 8125, 8215, 5145, 8725,
                                                 2887, 2977, 5845, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9400, 0, 3, 8215, 8305, 5245, 8875,
                                                 2977, 3067, 5995, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9625, 0, 3, 8305, 8395, 5345, 9025,
                                                 3067, 3157, 6145, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 9850, 0, 3, 8575, 8725, 5845, 9400,
                                                 3337, 3463, 6505, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 10165, 0, 3, 8725, 8875, 5995, 9625,
                                                 3463, 3589, 6715, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 10480, 0, 3, 9175, 9400, 6505, 10165,
                                                 3841, 4009, 7485, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10900, 3, 4345, 4355, 7780, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10921, 3, 4355, 4365, 7795, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10942, 3, 4365, 4375, 7810, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10963, 3, 4375, 4385, 7825, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10984, 3, 4385, 4395, 7840, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 11005, 0, 3, 7765, 10900, 7900, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11068, 0, 3, 7780, 10921, 7945, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11131, 0, 3, 7795, 10942, 7990, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11194, 0, 3, 7810, 10963, 8035, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11257, 0, 3, 7825, 10984, 8080, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 11320, 0, 3, 7855, 11005, 4625, 4685,
                                                 8125, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11446, 0, 3, 7900, 11068, 4685, 4745,
                                                 8215, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11572, 0, 3, 7945, 11131, 4745, 4805,
                                                 8305, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11698, 0, 3, 7990, 11194, 4805, 4865,
                                                 8395, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11824, 0, 3, 8035, 11257, 4865, 4925,
                                                 8485, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 11950, 0, 3, 8215, 11572, 5045, 5145,
                                                 8725, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12160, 0, 3, 8305, 11698, 5145, 5245,
                                                 8875, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12370, 0, 3, 8395, 11824, 5245, 5345,
                                                 9025, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 12580, 0, 3, 11320, 11446, 8575, 11950,
                                                 5545, 5695, 9175, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 12895, 0, 3, 11446, 11572, 8725, 12160,
                                                 5695, 5845, 9400, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 13210, 0, 3, 11572, 11698, 8875, 12370,
                                                 5845, 5995, 9625, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 13525, 0, 3, 11950, 12160, 9400, 13210,
                                                 6295, 6505, 10165, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 13966, 0, 3, 12580, 12895, 9850, 13525,
                                                 6925, 7205, 10480, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 14554, 13966, 588, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 15142, 14554, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 15142, 11, nmax);
}

}  // namespace simdt2ceri
