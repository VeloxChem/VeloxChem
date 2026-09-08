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


#include "SimdElectronRepulsionRecHH.hpp"

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

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 9298, nvalues);

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
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 17, 20, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 20, 23, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 23, 26, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 26, 29, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 29, 32, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 32, 35, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 35, 38, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 162, 0, 44, 50, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 177, 0, 50, 56, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 192, 0, 56, 62, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 207, 0, 62, 68, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 222, 0, 68, 74, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 237, 0, 74, 80, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 252, 0, 92, 102, 177, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 273, 0, 102, 112, 192, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 294, 0, 112, 122, 207, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 315, 0, 122, 132, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 336, 0, 132, 142, 237, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 357, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 360, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 363, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 366, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 369, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 372, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 375, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 378, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 387, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 396, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 405, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 414, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 423, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 432, 3, 15, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 441, 0, 3, 20, 378, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 459, 0, 3, 23, 387, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 477, 0, 3, 26, 396, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 495, 0, 3, 29, 405, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 513, 0, 3, 32, 414, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 531, 0, 3, 35, 423, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 549, 0, 3, 38, 432, 86, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 567, 0, 3, 44, 441, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 597, 0, 3, 50, 459, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 627, 0, 3, 56, 477, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 657, 0, 3, 62, 495, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 687, 0, 3, 68, 513, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 717, 0, 3, 74, 531, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 747, 0, 3, 80, 549, 152, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 777, 0, 3, 102, 627, 177, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 822, 0, 3, 112, 657, 192, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 867, 0, 3, 122, 687, 207, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 912, 0, 3, 132, 717, 222, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 957, 0, 3, 142, 747, 237, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1002, 0, 3, 162, 777, 252, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1065, 0, 3, 177, 822, 273, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1128, 0, 3, 192, 867, 294, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1191, 0, 3, 207, 912, 315, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1254, 0, 3, 222, 957, 336, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1317, 3, 9, 10, 360, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1323, 3, 10, 11, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1329, 3, 11, 12, 366, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1335, 3, 12, 13, 369, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1341, 3, 13, 14, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1347, 3, 14, 15, 375, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1353, 0, 3, 357, 1317, 387, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1371, 0, 3, 360, 1323, 396, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1389, 0, 3, 363, 1329, 405, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1407, 0, 3, 366, 1335, 414, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1425, 0, 3, 369, 1341, 423, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1443, 0, 3, 372, 1347, 432, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1461, 0, 3, 378, 1353, 44, 50, 459,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1497, 0, 3, 387, 1371, 50, 56, 477,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1533, 0, 3, 396, 1389, 56, 62, 495,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1569, 0, 3, 405, 1407, 62, 68, 513,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1605, 0, 3, 414, 1425, 68, 74, 531,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1641, 0, 3, 423, 1443, 74, 80, 549,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1677, 0, 3, 459, 1497, 92, 102, 627,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1737, 0, 3, 477, 1533, 102, 112, 657,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1797, 0, 3, 495, 1569, 112, 122, 687,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1857, 0, 3, 513, 1605, 122, 132, 717,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1917, 0, 3, 531, 1641, 132, 142, 747,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1977, 0, 3, 1461, 1497, 627, 1737, 162,
                                                 177, 822, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2067, 0, 3, 1497, 1533, 657, 1797, 177,
                                                 192, 867, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2157, 0, 3, 1533, 1569, 687, 1857, 192,
                                                 207, 912, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2247, 0, 3, 1569, 1605, 717, 1917, 207,
                                                 222, 957, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2337, 0, 3, 1677, 1737, 822, 2067, 252,
                                                 273, 1128, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2463, 0, 3, 1737, 1797, 867, 2157, 273,
                                                 294, 1191, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2589, 0, 3, 1797, 1857, 912, 2247, 294,
                                                 315, 1254, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2715, 3, 357, 360, 1323, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2725, 3, 360, 363, 1329, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2735, 3, 363, 366, 1335, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2745, 3, 366, 369, 1341, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2755, 3, 369, 372, 1347, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2765, 0, 3, 1317, 2715, 1371, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2795, 0, 3, 1323, 2725, 1389, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2825, 0, 3, 1329, 2735, 1407, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2855, 0, 3, 1335, 2745, 1425, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2885, 0, 3, 1341, 2755, 1443, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2915, 0, 3, 1353, 2765, 441, 459, 1497,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2975, 0, 3, 1371, 2795, 459, 477, 1533,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3035, 0, 3, 1389, 2825, 477, 495, 1569,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3095, 0, 3, 1407, 2855, 495, 513, 1605,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3155, 0, 3, 1425, 2885, 513, 531, 1641,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3215, 0, 3, 1461, 2915, 567, 597, 1677,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3315, 0, 3, 1497, 2975, 597, 627, 1737,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3415, 0, 3, 1533, 3035, 627, 657, 1797,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3515, 0, 3, 1569, 3095, 657, 687, 1857,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3615, 0, 3, 1605, 3155, 687, 717, 1917,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3715, 0, 3, 2915, 2975, 1737, 3415, 777,
                                                 822, 2067, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3865, 0, 3, 2975, 3035, 1797, 3515, 822,
                                                 867, 2157, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4015, 0, 3, 3035, 3095, 1857, 3615, 867,
                                                 912, 2247, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 4165, 0, 3, 3215, 3315, 1977, 3715,
                                                 1002, 1065, 2337, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 4375, 0, 3, 3315, 3415, 2067, 3865,
                                                 1065, 1128, 2463, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 4585, 0, 3, 3415, 3515, 2157, 4015,
                                                 1128, 1191, 2589, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4795, 3, 1317, 1323, 2725, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4810, 3, 1323, 1329, 2735, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4825, 3, 1329, 1335, 2745, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4840, 3, 1335, 1341, 2755, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4855, 0, 3, 2715, 4795, 2795, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4900, 0, 3, 2725, 4810, 2825, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4945, 0, 3, 2735, 4825, 2855, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4990, 0, 3, 2745, 4840, 2885, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 5035, 0, 3, 2765, 4855, 1461, 1497,
                                                 2975, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5125, 0, 3, 2795, 4900, 1497, 1533,
                                                 3035, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5215, 0, 3, 2825, 4945, 1533, 1569,
                                                 3095, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5305, 0, 3, 2855, 4990, 1569, 1605,
                                                 3155, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5395, 0, 3, 2975, 5125, 1677, 1737,
                                                 3415, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5545, 0, 3, 3035, 5215, 1737, 1797,
                                                 3515, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5695, 0, 3, 3095, 5305, 1797, 1857,
                                                 3615, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5845, 0, 3, 5035, 5125, 3415, 5545,
                                                 1977, 2067, 3865, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 6070, 0, 3, 5125, 5215, 3515, 5695,
                                                 2067, 2157, 4015, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 6295, 0, 3, 5395, 5545, 3865, 6070,
                                                 2337, 2463, 4585, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6610, 3, 2715, 2725, 4810, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6631, 3, 2725, 2735, 4825, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6652, 3, 2735, 2745, 4840, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 6673, 0, 3, 4795, 6610, 4900, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6736, 0, 3, 4810, 6631, 4945, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6799, 0, 3, 4825, 6652, 4990, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 6862, 0, 3, 4855, 6673, 2915, 2975,
                                                 5125, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6988, 0, 3, 4900, 6736, 2975, 3035,
                                                 5215, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7114, 0, 3, 4945, 6799, 3035, 3095,
                                                 5305, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7240, 0, 3, 5035, 6862, 3215, 3315,
                                                 5395, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7450, 0, 3, 5125, 6988, 3315, 3415,
                                                 5545, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7660, 0, 3, 5215, 7114, 3415, 3515,
                                                 5695, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 7870, 0, 3, 6862, 6988, 5545, 7660,
                                                 3715, 3865, 6070, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 8185, 0, 3, 7240, 7450, 5845, 7870,
                                                 4165, 4375, 6295, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8626, 8185, 441, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 9067, 8626, 21, nmax);

    simdtrf::transform_h_outer_tri(values, nvalues, buffer, 9067, nmax);
}

}  // namespace simdt2ceri
