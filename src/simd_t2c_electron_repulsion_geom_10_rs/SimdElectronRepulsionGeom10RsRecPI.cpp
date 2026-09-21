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


#include "SimdElectronRepulsionGeom10RsRecPI.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5059, 4516, 504, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 152, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 155, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 179, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 182, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 185, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 188, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 191, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 194, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 203, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 212, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 221, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 230, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 239, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 248, 3, 19, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 257, 3, 20, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 266, 3, 21, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 275, 3, 22, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 284, 3, 23, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 293, 3, 24, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 302, 0, 3, 29, 203, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 320, 0, 3, 32, 212, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 338, 0, 3, 35, 221, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 356, 0, 3, 38, 230, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 374, 0, 3, 41, 239, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 392, 0, 3, 50, 257, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 410, 0, 3, 53, 266, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 428, 0, 3, 56, 275, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 446, 0, 3, 59, 284, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 464, 0, 3, 62, 293, 146, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 482, 3, 7, 8, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 488, 3, 8, 9, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 494, 3, 9, 10, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 500, 3, 10, 11, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 506, 3, 11, 12, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 512, 3, 12, 13, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 518, 3, 13, 14, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 524, 3, 17, 18, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 530, 3, 18, 19, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 536, 3, 19, 20, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 542, 3, 20, 21, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 548, 3, 21, 22, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 554, 3, 22, 23, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 560, 3, 23, 24, 191, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 566, 0, 3, 155, 494, 203, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 584, 0, 3, 158, 500, 212, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 602, 0, 3, 161, 506, 221, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 620, 0, 3, 164, 512, 230, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 638, 0, 3, 167, 518, 239, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 656, 0, 3, 176, 536, 257, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 674, 0, 3, 179, 542, 266, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 692, 0, 3, 182, 548, 275, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 710, 0, 3, 185, 554, 284, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 728, 0, 3, 188, 560, 293, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 746, 0, 3, 194, 566, 68, 74, 302, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 782, 0, 3, 203, 584, 74, 80, 320, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 818, 0, 3, 212, 602, 80, 86, 338, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 854, 0, 3, 221, 620, 86, 92, 356, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 890, 0, 3, 230, 638, 92, 98, 374, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 926, 0, 3, 248, 656, 110, 116, 392,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 962, 0, 3, 257, 674, 116, 122, 410,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 998, 0, 3, 266, 692, 122, 128, 428,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1034, 0, 3, 275, 710, 128, 134, 446,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1070, 0, 3, 284, 728, 134, 140, 464,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1106, 3, 152, 155, 494, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1116, 3, 155, 158, 500, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1126, 3, 158, 161, 506, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1136, 3, 161, 164, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1146, 3, 164, 167, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1156, 3, 173, 176, 536, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1166, 3, 176, 179, 542, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1176, 3, 179, 182, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1186, 3, 182, 185, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1196, 3, 185, 188, 560, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1206, 0, 3, 494, 1116, 584, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1236, 0, 3, 500, 1126, 602, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1266, 0, 3, 506, 1136, 620, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1296, 0, 3, 512, 1146, 638, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1326, 0, 3, 536, 1166, 674, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1356, 0, 3, 542, 1176, 692, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1386, 0, 3, 548, 1186, 710, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1416, 0, 3, 554, 1196, 728, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1446, 0, 3, 584, 1236, 302, 320, 818,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1506, 0, 3, 602, 1266, 320, 338, 854,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1566, 0, 3, 620, 1296, 338, 356, 890,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1626, 0, 3, 674, 1356, 392, 410, 998,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1686, 0, 3, 692, 1386, 410, 428, 1034,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1746, 0, 3, 710, 1416, 428, 446, 1070,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1806, 3, 482, 488, 1106, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1821, 3, 488, 494, 1116, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1836, 3, 494, 500, 1126, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1851, 3, 500, 506, 1136, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1866, 3, 506, 512, 1146, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1881, 3, 524, 530, 1156, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1896, 3, 530, 536, 1166, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1911, 3, 536, 542, 1176, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1926, 3, 542, 548, 1186, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1941, 3, 548, 554, 1196, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1956, 0, 3, 1116, 1836, 1236, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2001, 0, 3, 1126, 1851, 1266, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2046, 0, 3, 1136, 1866, 1296, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2091, 0, 3, 1166, 1911, 1356, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2136, 0, 3, 1176, 1926, 1386, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2181, 0, 3, 1186, 1941, 1416, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2226, 0, 3, 1206, 1956, 746, 782, 1446,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2316, 0, 3, 1236, 2001, 782, 818, 1506,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2406, 0, 3, 1266, 2046, 818, 854, 1566,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2496, 0, 3, 1326, 2091, 926, 962, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2586, 0, 3, 1356, 2136, 962, 998, 1686,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2676, 0, 3, 1386, 2181, 998, 1034, 1746,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2766, 3, 1106, 1116, 1836, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2787, 3, 1116, 1126, 1851, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2808, 3, 1126, 1136, 1866, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2829, 3, 1156, 1166, 1911, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2850, 3, 1166, 1176, 1926, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2871, 3, 1176, 1186, 1941, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2892, 0, 3, 1836, 2787, 2001, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2955, 0, 3, 1851, 2808, 2046, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3018, 0, 3, 1911, 2850, 2136, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3081, 0, 3, 1926, 2871, 2181, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3144, 0, 3, 2001, 2955, 1446, 1506,
                                                 2406, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3270, 0, 3, 2136, 3081, 1626, 1686,
                                                 2676, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3396, 3, 1806, 1821, 2766, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3424, 3, 1836, 1851, 2808, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3452, 3, 1881, 1896, 2829, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3480, 3, 1911, 1926, 2871, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 3508, 0, 3, 2787, 3424, 2955, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 3592, 0, 3, 2850, 3480, 3081, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3676, 0, 3, 2892, 3508, 2226, 2316,
                                                 3144, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3844, 0, 3, 3018, 3592, 2496, 2586,
                                                 3270, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 4012, 3452, 3844, 1, 28, ncols, alpha);

            simdgeo::geom_p_y(buffer, 4096, 3452, 3844, 1, 28, ncols, alpha);

            simdgeo::geom_p_z(buffer, 4180, 3452, 3844, 1, 28, ncols, alpha);

            simdgeo::geom_p_x(buffer, 4264, 3396, 3676, 1, 28, ncols, alpha);

            simdgeo::geom_p_y(buffer, 4348, 3396, 3676, 1, 28, ncols, alpha);

            simdgeo::geom_p_z(buffer, 4432, 3396, 3676, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 4516, 4264, 252, ncols);

            simdfunc::contract_primitives(buffer, 4768, 4012, 252, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 5020, 4768, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 5020, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5020, 4852, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 39 * nvalues, nvalues, buffer, 5020, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5020, 4936, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 78 * nvalues, nvalues, buffer, 5020, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5020, 4516, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 117 * nvalues, nvalues, buffer, 5020, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5020, 4600, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 156 * nvalues, nvalues, buffer, 5020, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5020, 4684, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 195 * nvalues, nvalues, buffer, 5020, 13, nmax);
}

}  // namespace simdt2ceri
