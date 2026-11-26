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


#include "SimdElectronRepulsionGeom10RsRecPK.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 7657, 6964, 648, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 179, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 182, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 185, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 188, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 191, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 194, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 197, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 200, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 203, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 206, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 209, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 212, 3, 8, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 221, 3, 9, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 230, 3, 10, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 239, 3, 11, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 248, 3, 12, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 257, 3, 13, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 266, 3, 14, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 275, 3, 18, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 284, 3, 19, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 293, 3, 20, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 302, 3, 21, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 311, 3, 22, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 320, 3, 23, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 329, 3, 24, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 338, 0, 3, 26, 212, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 356, 0, 3, 29, 221, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 374, 0, 3, 32, 230, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 392, 0, 3, 35, 239, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 410, 0, 3, 38, 248, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 428, 0, 3, 41, 257, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 446, 0, 3, 44, 266, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 464, 0, 3, 50, 275, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 482, 0, 3, 53, 284, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 500, 0, 3, 56, 293, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 518, 0, 3, 59, 302, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 536, 0, 3, 62, 311, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 554, 0, 3, 65, 320, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 572, 0, 3, 68, 329, 152, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 590, 3, 7, 8, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 596, 3, 8, 9, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 602, 3, 9, 10, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 608, 3, 10, 11, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 614, 3, 11, 12, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 620, 3, 12, 13, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 626, 3, 13, 14, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 632, 3, 17, 18, 191, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 638, 3, 18, 19, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 644, 3, 19, 20, 197, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 650, 3, 20, 21, 200, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 656, 3, 21, 22, 203, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 662, 3, 22, 23, 206, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 668, 3, 23, 24, 209, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 674, 0, 3, 164, 596, 221, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 692, 0, 3, 167, 602, 230, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 710, 0, 3, 170, 608, 239, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 728, 0, 3, 173, 614, 248, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 746, 0, 3, 176, 620, 257, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 764, 0, 3, 179, 626, 266, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 782, 0, 3, 191, 638, 284, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 800, 0, 3, 194, 644, 293, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 818, 0, 3, 197, 650, 302, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 836, 0, 3, 200, 656, 311, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 854, 0, 3, 203, 662, 320, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 872, 0, 3, 206, 668, 329, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 890, 0, 3, 221, 692, 74, 80, 374, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 926, 0, 3, 230, 710, 80, 86, 392, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 962, 0, 3, 239, 728, 86, 92, 410, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 998, 0, 3, 248, 746, 92, 98, 428, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1034, 0, 3, 257, 764, 98, 104, 446,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1070, 0, 3, 284, 800, 116, 122, 500,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1106, 0, 3, 293, 818, 122, 128, 518,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1142, 0, 3, 302, 836, 128, 134, 536,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1178, 0, 3, 311, 854, 134, 140, 554,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1214, 0, 3, 320, 872, 140, 146, 572,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1250, 3, 158, 161, 590, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1260, 3, 161, 164, 596, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1270, 3, 164, 167, 602, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1280, 3, 167, 170, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1290, 3, 170, 173, 614, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1300, 3, 173, 176, 620, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1310, 3, 176, 179, 626, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1320, 3, 185, 188, 632, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1330, 3, 188, 191, 638, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1340, 3, 191, 194, 644, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1350, 3, 194, 197, 650, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1360, 3, 197, 200, 656, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1370, 3, 200, 203, 662, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1380, 3, 203, 206, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1390, 0, 3, 596, 1270, 692, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1420, 0, 3, 602, 1280, 710, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1450, 0, 3, 608, 1290, 728, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1480, 0, 3, 614, 1300, 746, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1510, 0, 3, 620, 1310, 764, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1540, 0, 3, 638, 1340, 800, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1570, 0, 3, 644, 1350, 818, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1600, 0, 3, 650, 1360, 836, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1630, 0, 3, 656, 1370, 854, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1660, 0, 3, 662, 1380, 872, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1690, 0, 3, 674, 1390, 338, 356, 890,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1750, 0, 3, 692, 1420, 356, 374, 926,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1810, 0, 3, 710, 1450, 374, 392, 962,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1870, 0, 3, 728, 1480, 392, 410, 998,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1930, 0, 3, 746, 1510, 410, 428, 1034,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1990, 0, 3, 782, 1540, 464, 482, 1070,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2050, 0, 3, 800, 1570, 482, 500, 1106,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2110, 0, 3, 818, 1600, 500, 518, 1142,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2170, 0, 3, 836, 1630, 518, 536, 1178,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2230, 0, 3, 854, 1660, 536, 554, 1214,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2290, 3, 590, 596, 1270, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2305, 3, 596, 602, 1280, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2320, 3, 602, 608, 1290, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2335, 3, 608, 614, 1300, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2350, 3, 614, 620, 1310, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2365, 3, 632, 638, 1340, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2380, 3, 638, 644, 1350, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2395, 3, 644, 650, 1360, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2410, 3, 650, 656, 1370, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2425, 3, 656, 662, 1380, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2440, 0, 3, 1270, 2305, 1420, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2485, 0, 3, 1280, 2320, 1450, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2530, 0, 3, 1290, 2335, 1480, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2575, 0, 3, 1300, 2350, 1510, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2620, 0, 3, 1340, 2380, 1570, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2665, 0, 3, 1350, 2395, 1600, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2710, 0, 3, 1360, 2410, 1630, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2755, 0, 3, 1370, 2425, 1660, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2800, 0, 3, 1420, 2485, 890, 926, 1810,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2890, 0, 3, 1450, 2530, 926, 962, 1870,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2980, 0, 3, 1480, 2575, 962, 998, 1930,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3070, 0, 3, 1570, 2665, 1070, 1106,
                                                 2110, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3160, 0, 3, 1600, 2710, 1106, 1142,
                                                 2170, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3250, 0, 3, 1630, 2755, 1142, 1178,
                                                 2230, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3340, 3, 1250, 1260, 2290, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3361, 3, 1260, 1270, 2305, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3382, 3, 1270, 1280, 2320, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3403, 3, 1280, 1290, 2335, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3424, 3, 1290, 1300, 2350, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3445, 3, 1320, 1330, 2365, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3466, 3, 1330, 1340, 2380, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3487, 3, 1340, 1350, 2395, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3508, 3, 1350, 1360, 2410, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3529, 3, 1360, 1370, 2425, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3550, 0, 3, 2305, 3382, 2485, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3613, 0, 3, 2320, 3403, 2530, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3676, 0, 3, 2335, 3424, 2575, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3739, 0, 3, 2380, 3487, 2665, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3802, 0, 3, 2395, 3508, 2710, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3865, 0, 3, 2410, 3529, 2755, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3928, 0, 3, 2440, 3550, 1690, 1750,
                                                 2800, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4054, 0, 3, 2485, 3613, 1750, 1810,
                                                 2890, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4180, 0, 3, 2530, 3676, 1810, 1870,
                                                 2980, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4306, 0, 3, 2620, 3739, 1990, 2050,
                                                 3070, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4432, 0, 3, 2665, 3802, 2050, 2110,
                                                 3160, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4558, 0, 3, 2710, 3865, 2110, 2170,
                                                 3250, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4684, 3, 2290, 2305, 3382, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4712, 3, 2305, 2320, 3403, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4740, 3, 2320, 2335, 3424, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4768, 3, 2365, 2380, 3487, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4796, 3, 2380, 2395, 3508, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4824, 3, 2395, 2410, 3529, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 4852, 0, 3, 3382, 4712, 3613, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4936, 0, 3, 3403, 4740, 3676, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5020, 0, 3, 3487, 4796, 3802, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5104, 0, 3, 3508, 4824, 3865, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 5188, 0, 3, 3613, 4936, 2800, 2890,
                                                 4180, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 5356, 0, 3, 3802, 5104, 3070, 3160,
                                                 4558, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5524, 3, 3340, 3361, 4684, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5560, 3, 3382, 3403, 4740, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5596, 3, 3445, 3466, 4768, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5632, 3, 3487, 3508, 4824, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 5668, 0, 3, 4712, 5560, 4936, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 5776, 0, 3, 4796, 5632, 5104, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 5884, 0, 3, 4852, 5668, 3928, 4054,
                                                 5188, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 6100, 0, 3, 5020, 5776, 4306, 4432,
                                                 5356, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 6316, 5596, 6100, 1, 36, ncols, alpha);

            simdgeo::geom_p_y(buffer, 6424, 5596, 6100, 1, 36, ncols, alpha);

            simdgeo::geom_p_z(buffer, 6532, 5596, 6100, 1, 36, ncols, alpha);

            simdgeo::geom_p_x(buffer, 6640, 5524, 5884, 1, 36, ncols, alpha);

            simdgeo::geom_p_y(buffer, 6748, 5524, 5884, 1, 36, ncols, alpha);

            simdgeo::geom_p_z(buffer, 6856, 5524, 5884, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 6964, 6640, 324, ncols);

            simdfunc::contract_primitives(buffer, 7288, 6316, 324, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 7612, 7288, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 7612, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7612, 7396, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 45 * nvalues, nvalues, buffer, 7612, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7612, 7504, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 90 * nvalues, nvalues, buffer, 7612, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7612, 6964, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 135 * nvalues, nvalues, buffer, 7612, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7612, 7072, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 180 * nvalues, nvalues, buffer, 7612, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7612, 7180, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 225 * nvalues, nvalues, buffer, 7612, 15, nmax);
}

}  // namespace simdt2ceri
