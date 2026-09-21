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


#include "SimdElectronRepulsionRsRecDI.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_di_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_di_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4218, 3804, 336, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 8, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 16, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 9, 10, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 10, 11, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 11, 12, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 12, 13, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 13, 14, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 17, 18, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 18, 19, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 19, 20, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 20, 21, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 21, 22, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 22, 23, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 23, 24, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 152, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 155, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 179, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 182, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 185, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 188, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 197, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 206, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 215, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 224, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 233, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 242, 3, 19, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 251, 3, 20, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 260, 3, 21, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 269, 3, 22, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 278, 3, 23, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 287, 3, 24, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 296, 0, 3, 29, 197, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 314, 0, 3, 32, 206, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 332, 0, 3, 35, 215, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 350, 0, 3, 38, 224, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 368, 0, 3, 41, 233, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 386, 0, 3, 50, 251, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 404, 0, 3, 53, 260, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 422, 0, 3, 56, 269, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 440, 0, 3, 59, 278, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 458, 0, 3, 62, 287, 146, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 476, 3, 9, 10, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 482, 3, 10, 11, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 488, 3, 11, 12, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 494, 3, 12, 13, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 500, 3, 13, 14, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 506, 3, 19, 20, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 512, 3, 20, 21, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 518, 3, 21, 22, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 524, 3, 22, 23, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 530, 3, 23, 24, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 536, 0, 3, 152, 476, 197, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 554, 0, 3, 155, 482, 206, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 572, 0, 3, 158, 488, 215, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 590, 0, 3, 161, 494, 224, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 608, 0, 3, 164, 500, 233, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 626, 0, 3, 170, 506, 251, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 644, 0, 3, 173, 512, 260, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 662, 0, 3, 176, 518, 269, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 680, 0, 3, 179, 524, 278, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 698, 0, 3, 182, 530, 287, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 716, 0, 3, 188, 536, 68, 74, 296, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 752, 0, 3, 197, 554, 74, 80, 314, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 788, 0, 3, 206, 572, 80, 86, 332, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 824, 0, 3, 215, 590, 86, 92, 350, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 860, 0, 3, 224, 608, 92, 98, 368, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 896, 0, 3, 242, 626, 110, 116, 386,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 932, 0, 3, 251, 644, 116, 122, 404,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 968, 0, 3, 260, 662, 122, 128, 422,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1004, 0, 3, 269, 680, 128, 134, 440,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1040, 0, 3, 278, 698, 134, 140, 458,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1076, 3, 152, 155, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1086, 3, 155, 158, 488, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1096, 3, 158, 161, 494, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1106, 3, 161, 164, 500, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1116, 3, 170, 173, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1126, 3, 173, 176, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1136, 3, 176, 179, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1146, 3, 179, 182, 530, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1156, 0, 3, 476, 1076, 554, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1186, 0, 3, 482, 1086, 572, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1216, 0, 3, 488, 1096, 590, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1246, 0, 3, 494, 1106, 608, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1276, 0, 3, 506, 1116, 644, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1306, 0, 3, 512, 1126, 662, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1336, 0, 3, 518, 1136, 680, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1366, 0, 3, 524, 1146, 698, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1396, 0, 3, 554, 1186, 296, 314, 788,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1456, 0, 3, 572, 1216, 314, 332, 824,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1516, 0, 3, 590, 1246, 332, 350, 860,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1576, 0, 3, 644, 1306, 386, 404, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1636, 0, 3, 662, 1336, 404, 422, 1004,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1696, 0, 3, 680, 1366, 422, 440, 1040,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1756, 3, 476, 482, 1086, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1771, 3, 482, 488, 1096, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1786, 3, 488, 494, 1106, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1801, 3, 506, 512, 1126, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1816, 3, 512, 518, 1136, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1831, 3, 518, 524, 1146, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1846, 0, 3, 1076, 1756, 1186, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1891, 0, 3, 1086, 1771, 1216, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1936, 0, 3, 1096, 1786, 1246, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1981, 0, 3, 1116, 1801, 1306, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2026, 0, 3, 1126, 1816, 1336, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2071, 0, 3, 1136, 1831, 1366, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2116, 0, 3, 1156, 1846, 716, 752, 1396,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2206, 0, 3, 1186, 1891, 752, 788, 1456,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2296, 0, 3, 1216, 1936, 788, 824, 1516,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2386, 0, 3, 1276, 1981, 896, 932, 1576,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2476, 0, 3, 1306, 2026, 932, 968, 1636,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2566, 0, 3, 1336, 2071, 968, 1004, 1696,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2656, 3, 1076, 1086, 1771, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2677, 3, 1086, 1096, 1786, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2698, 3, 1116, 1126, 1816, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2719, 3, 1126, 1136, 1831, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2740, 0, 3, 1756, 2656, 1891, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2803, 0, 3, 1771, 2677, 1936, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2866, 0, 3, 1801, 2698, 2026, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2929, 0, 3, 1816, 2719, 2071, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2992, 0, 3, 1891, 2803, 1396, 1456,
                                                 2296, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3118, 0, 3, 2026, 2929, 1576, 1636,
                                                 2566, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3244, 3, 1756, 1771, 2677, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3272, 3, 1801, 1816, 2719, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 3300, 0, 3, 2656, 3244, 2803, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 3384, 0, 3, 2698, 3272, 2929, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3468, 0, 3, 2740, 3300, 2116, 2206,
                                                 2992, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3636, 0, 3, 2866, 3384, 2386, 2476,
                                                 3118, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3804, 3468, 336, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 4140, 3972, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 4140, 13, nmax);

    simdtrf::transform_i_inner(buffer, 4140, 3804, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 65 * nvalues, nvalues, buffer, 4140, 13, nmax);
}

}  // namespace simdt2ceri
