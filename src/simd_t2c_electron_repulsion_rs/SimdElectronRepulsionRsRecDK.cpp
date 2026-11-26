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


#include "SimdElectronRepulsionRsRecDK.hpp"

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
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_dk_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_dk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 6532, 6010, 432, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 16, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 7, 8, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 8, 9, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 9, 10, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 10, 11, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 11, 12, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 12, 13, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 13, 14, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 18, 19, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 19, 20, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 20, 21, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 21, 22, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 22, 23, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 23, 24, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 179, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 182, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 185, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 188, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 191, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 194, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 197, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 200, 3, 8, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 209, 3, 9, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 218, 3, 10, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 227, 3, 11, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 236, 3, 12, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 245, 3, 13, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 254, 3, 14, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 263, 3, 18, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 272, 3, 19, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 281, 3, 20, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 290, 3, 21, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 299, 3, 22, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 308, 3, 23, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 317, 3, 24, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 326, 0, 3, 26, 200, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 344, 0, 3, 29, 209, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 362, 0, 3, 32, 218, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 380, 0, 3, 35, 227, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 398, 0, 3, 38, 236, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 416, 0, 3, 41, 245, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 434, 0, 3, 44, 254, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 452, 0, 3, 50, 263, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 470, 0, 3, 53, 272, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 488, 0, 3, 56, 281, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 506, 0, 3, 59, 290, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 524, 0, 3, 62, 299, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 542, 0, 3, 65, 308, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 560, 0, 3, 68, 317, 152, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 578, 3, 8, 9, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 584, 3, 9, 10, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 590, 3, 10, 11, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 596, 3, 11, 12, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 602, 3, 12, 13, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 608, 3, 13, 14, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 614, 3, 18, 19, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 620, 3, 19, 20, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 626, 3, 20, 21, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 632, 3, 21, 22, 191, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 638, 3, 22, 23, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 644, 3, 23, 24, 197, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 650, 0, 3, 158, 578, 209, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 668, 0, 3, 161, 584, 218, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 686, 0, 3, 164, 590, 227, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 704, 0, 3, 167, 596, 236, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 722, 0, 3, 170, 602, 245, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 740, 0, 3, 173, 608, 254, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 758, 0, 3, 179, 614, 272, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 776, 0, 3, 182, 620, 281, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 794, 0, 3, 185, 626, 290, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 812, 0, 3, 188, 632, 299, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 830, 0, 3, 191, 638, 308, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 848, 0, 3, 194, 644, 317, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 866, 0, 3, 209, 668, 74, 80, 362, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 902, 0, 3, 218, 686, 80, 86, 380, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 938, 0, 3, 227, 704, 86, 92, 398, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 974, 0, 3, 236, 722, 92, 98, 416, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1010, 0, 3, 245, 740, 98, 104, 434,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1046, 0, 3, 272, 776, 116, 122, 488,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1082, 0, 3, 281, 794, 122, 128, 506,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1118, 0, 3, 290, 812, 128, 134, 524,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1154, 0, 3, 299, 830, 134, 140, 542,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1190, 0, 3, 308, 848, 140, 146, 560,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1226, 3, 158, 161, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1236, 3, 161, 164, 590, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1246, 3, 164, 167, 596, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1256, 3, 167, 170, 602, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1266, 3, 170, 173, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1276, 3, 179, 182, 620, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1286, 3, 182, 185, 626, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1296, 3, 185, 188, 632, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1306, 3, 188, 191, 638, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1316, 3, 191, 194, 644, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1326, 0, 3, 578, 1226, 668, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1356, 0, 3, 584, 1236, 686, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1386, 0, 3, 590, 1246, 704, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1416, 0, 3, 596, 1256, 722, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1446, 0, 3, 602, 1266, 740, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1476, 0, 3, 614, 1276, 776, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1506, 0, 3, 620, 1286, 794, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1536, 0, 3, 626, 1296, 812, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1566, 0, 3, 632, 1306, 830, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1596, 0, 3, 638, 1316, 848, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1626, 0, 3, 650, 1326, 326, 344, 866,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1686, 0, 3, 668, 1356, 344, 362, 902,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1746, 0, 3, 686, 1386, 362, 380, 938,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1806, 0, 3, 704, 1416, 380, 398, 974,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1866, 0, 3, 722, 1446, 398, 416, 1010,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1926, 0, 3, 758, 1476, 452, 470, 1046,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1986, 0, 3, 776, 1506, 470, 488, 1082,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2046, 0, 3, 794, 1536, 488, 506, 1118,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2106, 0, 3, 812, 1566, 506, 524, 1154,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2166, 0, 3, 830, 1596, 524, 542, 1190,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2226, 3, 578, 584, 1236, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2241, 3, 584, 590, 1246, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2256, 3, 590, 596, 1256, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2271, 3, 596, 602, 1266, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2286, 3, 614, 620, 1286, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2301, 3, 620, 626, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2316, 3, 626, 632, 1306, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2331, 3, 632, 638, 1316, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2346, 0, 3, 1226, 2226, 1356, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2391, 0, 3, 1236, 2241, 1386, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2436, 0, 3, 1246, 2256, 1416, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2481, 0, 3, 1256, 2271, 1446, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2526, 0, 3, 1276, 2286, 1506, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2571, 0, 3, 1286, 2301, 1536, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2616, 0, 3, 1296, 2316, 1566, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2661, 0, 3, 1306, 2331, 1596, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2706, 0, 3, 1356, 2391, 866, 902, 1746,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2796, 0, 3, 1386, 2436, 902, 938, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2886, 0, 3, 1416, 2481, 938, 974, 1866,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2976, 0, 3, 1506, 2571, 1046, 1082,
                                                 2046, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3066, 0, 3, 1536, 2616, 1082, 1118,
                                                 2106, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3156, 0, 3, 1566, 2661, 1118, 1154,
                                                 2166, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3246, 3, 1226, 1236, 2241, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3267, 3, 1236, 1246, 2256, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3288, 3, 1246, 1256, 2271, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3309, 3, 1276, 1286, 2301, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3330, 3, 1286, 1296, 2316, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3351, 3, 1296, 1306, 2331, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3372, 0, 3, 2226, 3246, 2391, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3435, 0, 3, 2241, 3267, 2436, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3498, 0, 3, 2256, 3288, 2481, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3561, 0, 3, 2286, 3309, 2571, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3624, 0, 3, 2301, 3330, 2616, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3687, 0, 3, 2316, 3351, 2661, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3750, 0, 3, 2346, 3372, 1626, 1686,
                                                 2706, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3876, 0, 3, 2391, 3435, 1686, 1746,
                                                 2796, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4002, 0, 3, 2436, 3498, 1746, 1806,
                                                 2886, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4128, 0, 3, 2526, 3561, 1926, 1986,
                                                 2976, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4254, 0, 3, 2571, 3624, 1986, 2046,
                                                 3066, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4380, 0, 3, 2616, 3687, 2046, 2106,
                                                 3156, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4506, 3, 2226, 2241, 3267, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4534, 3, 2241, 2256, 3288, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4562, 3, 2286, 2301, 3330, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4590, 3, 2301, 2316, 3351, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 4618, 0, 3, 3246, 4506, 3435, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4702, 0, 3, 3267, 4534, 3498, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4786, 0, 3, 3309, 4562, 3624, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4870, 0, 3, 3330, 4590, 3687, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 4954, 0, 3, 3435, 4702, 2706, 2796,
                                                 4002, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 5122, 0, 3, 3624, 4870, 2976, 3066,
                                                 4380, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5290, 3, 3246, 3267, 4534, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5326, 3, 3309, 3330, 4590, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 5362, 0, 3, 4506, 5290, 4702, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 5470, 0, 3, 4562, 5326, 4870, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 5578, 0, 3, 4618, 5362, 3750, 3876,
                                                 4954, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 5794, 0, 3, 4786, 5470, 4128, 4254,
                                                 5122, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 6010, 5578, 432, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 6442, 6226, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 6442, 15, nmax);

    simdtrf::transform_k_inner(buffer, 6442, 6010, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 75 * nvalues, nvalues, buffer, 6442, 15, nmax);
}

}  // namespace simdt2ceri
