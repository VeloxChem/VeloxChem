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


#include "SimdElectronRepulsionGeom10RsRecDL.hpp"

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
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFL.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_dl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 22818, 21096, 1620, nvalues);

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
                                                9, 10, 11}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 18, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 7, 8, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 8, 9, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 9, 10, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 10, 11, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 11, 12, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 12, 13, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 13, 14, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 14, 15, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 15, 16, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 19, 20, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 20, 21, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 21, 22, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 22, 23, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 23, 24, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 24, 25, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 25, 26, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 26, 27, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 27, 28, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 51, 54, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 54, 57, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 63, 66, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 304, 0, 66, 69, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 314, 0, 69, 72, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 72, 75, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 75, 78, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 78, 81, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 81, 84, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 84, 87, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 87, 90, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 384, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 387, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 390, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 393, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 396, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 399, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 402, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 405, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 408, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 411, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 414, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 417, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 420, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 423, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 426, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 429, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 432, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 435, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 438, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 441, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 444, 3, 9, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 453, 3, 10, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 462, 3, 11, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 471, 3, 12, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 480, 3, 13, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 489, 3, 14, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 498, 3, 15, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 507, 3, 16, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 516, 3, 21, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 525, 3, 22, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 534, 3, 23, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 543, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 552, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 561, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 570, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 579, 3, 28, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 588, 0, 3, 36, 444, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 606, 0, 3, 39, 453, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 624, 0, 3, 42, 462, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 642, 0, 3, 45, 471, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 660, 0, 3, 48, 480, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 678, 0, 3, 51, 489, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 696, 0, 3, 54, 498, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 714, 0, 3, 57, 507, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 732, 0, 3, 69, 516, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 750, 0, 3, 72, 525, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 768, 0, 3, 75, 534, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 786, 0, 3, 78, 543, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 804, 0, 3, 81, 552, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 822, 0, 3, 84, 561, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 840, 0, 3, 87, 570, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 858, 0, 3, 90, 579, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 876, 0, 3, 102, 606, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 906, 0, 3, 108, 624, 234, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 936, 0, 3, 114, 642, 244, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 966, 0, 3, 120, 660, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 996, 0, 3, 126, 678, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1026, 0, 3, 132, 696, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1056, 0, 3, 138, 714, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1086, 0, 3, 156, 750, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1116, 0, 3, 162, 768, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1146, 0, 3, 168, 786, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1176, 0, 3, 174, 804, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1206, 0, 3, 180, 822, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1236, 0, 3, 186, 840, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1266, 0, 3, 192, 858, 374, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1296, 3, 7, 8, 387, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1302, 3, 8, 9, 390, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1308, 3, 9, 10, 393, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1314, 3, 10, 11, 396, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1320, 3, 11, 12, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1326, 3, 12, 13, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1332, 3, 13, 14, 405, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1338, 3, 14, 15, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1344, 3, 15, 16, 411, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1350, 3, 19, 20, 417, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1356, 3, 20, 21, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1362, 3, 21, 22, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1368, 3, 22, 23, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1374, 3, 23, 24, 429, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1380, 3, 24, 25, 432, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1386, 3, 25, 26, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1392, 3, 26, 27, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1398, 3, 27, 28, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1404, 0, 3, 390, 1308, 453, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1422, 0, 3, 393, 1314, 462, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1440, 0, 3, 396, 1320, 471, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1458, 0, 3, 399, 1326, 480, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1476, 0, 3, 402, 1332, 489, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1494, 0, 3, 405, 1338, 498, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1512, 0, 3, 408, 1344, 507, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1530, 0, 3, 420, 1362, 525, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1548, 0, 3, 423, 1368, 534, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1566, 0, 3, 426, 1374, 543, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1584, 0, 3, 429, 1380, 552, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1602, 0, 3, 432, 1386, 561, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1620, 0, 3, 435, 1392, 570, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1638, 0, 3, 438, 1398, 579, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1656, 0, 3, 444, 1404, 96, 102, 606,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1692, 0, 3, 453, 1422, 102, 108, 624,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1728, 0, 3, 462, 1440, 108, 114, 642,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1764, 0, 3, 471, 1458, 114, 120, 660,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1800, 0, 3, 480, 1476, 120, 126, 678,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1836, 0, 3, 489, 1494, 126, 132, 696,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1872, 0, 3, 498, 1512, 132, 138, 714,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1908, 0, 3, 516, 1530, 150, 156, 750,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1944, 0, 3, 525, 1548, 156, 162, 768,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1980, 0, 3, 534, 1566, 162, 168, 786,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2016, 0, 3, 543, 1584, 168, 174, 804,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2052, 0, 3, 552, 1602, 174, 180, 822,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2088, 0, 3, 561, 1620, 180, 186, 840,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2124, 0, 3, 570, 1638, 186, 192, 858,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2160, 0, 3, 588, 1656, 204, 214, 876,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2220, 0, 3, 606, 1692, 214, 224, 906,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2280, 0, 3, 624, 1728, 224, 234, 936,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2340, 0, 3, 642, 1764, 234, 244, 966,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2400, 0, 3, 660, 1800, 244, 254, 996,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2460, 0, 3, 678, 1836, 254, 264, 1026,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2520, 0, 3, 696, 1872, 264, 274, 1056,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2580, 0, 3, 732, 1908, 294, 304, 1086,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2640, 0, 3, 750, 1944, 304, 314, 1116,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2700, 0, 3, 768, 1980, 314, 324, 1146,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2760, 0, 3, 786, 2016, 324, 334, 1176,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2820, 0, 3, 804, 2052, 334, 344, 1206,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2880, 0, 3, 822, 2088, 344, 354, 1236,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2940, 0, 3, 840, 2124, 354, 364, 1266,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3000, 3, 384, 387, 1302, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3010, 3, 387, 390, 1308, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3020, 3, 390, 393, 1314, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3030, 3, 393, 396, 1320, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3040, 3, 396, 399, 1326, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3050, 3, 399, 402, 1332, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3060, 3, 402, 405, 1338, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3070, 3, 405, 408, 1344, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3080, 3, 414, 417, 1356, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3090, 3, 417, 420, 1362, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3100, 3, 420, 423, 1368, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3110, 3, 423, 426, 1374, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3120, 3, 426, 429, 1380, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3130, 3, 429, 432, 1386, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3140, 3, 432, 435, 1392, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3150, 3, 435, 438, 1398, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3160, 0, 3, 1308, 3020, 1422, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3190, 0, 3, 1314, 3030, 1440, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3220, 0, 3, 1320, 3040, 1458, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3250, 0, 3, 1326, 3050, 1476, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3280, 0, 3, 1332, 3060, 1494, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3310, 0, 3, 1338, 3070, 1512, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3340, 0, 3, 1362, 3100, 1548, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3370, 0, 3, 1368, 3110, 1566, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3400, 0, 3, 1374, 3120, 1584, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3430, 0, 3, 1380, 3130, 1602, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3460, 0, 3, 1386, 3140, 1620, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3490, 0, 3, 1392, 3150, 1638, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3520, 0, 3, 1404, 3160, 588, 606, 1692,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3580, 0, 3, 1422, 3190, 606, 624, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3640, 0, 3, 1440, 3220, 624, 642, 1764,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3700, 0, 3, 1458, 3250, 642, 660, 1800,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3760, 0, 3, 1476, 3280, 660, 678, 1836,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3820, 0, 3, 1494, 3310, 678, 696, 1872,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3880, 0, 3, 1530, 3340, 732, 750, 1944,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3940, 0, 3, 1548, 3370, 750, 768, 1980,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4000, 0, 3, 1566, 3400, 768, 786, 2016,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4060, 0, 3, 1584, 3430, 786, 804, 2052,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4120, 0, 3, 1602, 3460, 804, 822, 2088,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4180, 0, 3, 1620, 3490, 822, 840, 2124,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4240, 0, 3, 1692, 3580, 876, 906, 2280,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4340, 0, 3, 1728, 3640, 906, 936, 2340,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4440, 0, 3, 1764, 3700, 936, 966, 2400,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4540, 0, 3, 1800, 3760, 966, 996, 2460,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4640, 0, 3, 1836, 3820, 996, 1026, 2520,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4740, 0, 3, 1944, 3940, 1086, 1116,
                                                 2700, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4840, 0, 3, 1980, 4000, 1116, 1146,
                                                 2760, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4940, 0, 3, 2016, 4060, 1146, 1176,
                                                 2820, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5040, 0, 3, 2052, 4120, 1176, 1206,
                                                 2880, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5140, 0, 3, 2088, 4180, 1206, 1236,
                                                 2940, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5240, 3, 1296, 1302, 3010, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5255, 3, 1302, 1308, 3020, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5270, 3, 1308, 1314, 3030, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5285, 3, 1314, 1320, 3040, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5300, 3, 1320, 1326, 3050, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5315, 3, 1326, 1332, 3060, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5330, 3, 1332, 1338, 3070, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5345, 3, 1350, 1356, 3090, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5360, 3, 1356, 1362, 3100, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5375, 3, 1362, 1368, 3110, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5390, 3, 1368, 1374, 3120, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5405, 3, 1374, 1380, 3130, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5420, 3, 1380, 1386, 3140, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5435, 3, 1386, 1392, 3150, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 5450, 0, 3, 3020, 5270, 3190, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5495, 0, 3, 3030, 5285, 3220, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5540, 0, 3, 3040, 5300, 3250, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5585, 0, 3, 3050, 5315, 3280, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5630, 0, 3, 3060, 5330, 3310, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5675, 0, 3, 3100, 5375, 3370, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5720, 0, 3, 3110, 5390, 3400, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5765, 0, 3, 3120, 5405, 3430, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5810, 0, 3, 3130, 5420, 3460, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5855, 0, 3, 3140, 5435, 3490, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 5900, 0, 3, 3160, 5450, 1656, 1692,
                                                 3580, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5990, 0, 3, 3190, 5495, 1692, 1728,
                                                 3640, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6080, 0, 3, 3220, 5540, 1728, 1764,
                                                 3700, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6170, 0, 3, 3250, 5585, 1764, 1800,
                                                 3760, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6260, 0, 3, 3280, 5630, 1800, 1836,
                                                 3820, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6350, 0, 3, 3340, 5675, 1908, 1944,
                                                 3940, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6440, 0, 3, 3370, 5720, 1944, 1980,
                                                 4000, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6530, 0, 3, 3400, 5765, 1980, 2016,
                                                 4060, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6620, 0, 3, 3430, 5810, 2016, 2052,
                                                 4120, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6710, 0, 3, 3460, 5855, 2052, 2088,
                                                 4180, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6800, 0, 3, 3520, 5900, 2160, 2220,
                                                 4240, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6950, 0, 3, 3580, 5990, 2220, 2280,
                                                 4340, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7100, 0, 3, 3640, 6080, 2280, 2340,
                                                 4440, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7250, 0, 3, 3700, 6170, 2340, 2400,
                                                 4540, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7400, 0, 3, 3760, 6260, 2400, 2460,
                                                 4640, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7550, 0, 3, 3880, 6350, 2580, 2640,
                                                 4740, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7700, 0, 3, 3940, 6440, 2640, 2700,
                                                 4840, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7850, 0, 3, 4000, 6530, 2700, 2760,
                                                 4940, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8000, 0, 3, 4060, 6620, 2760, 2820,
                                                 5040, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8150, 0, 3, 4120, 6710, 2820, 2880,
                                                 5140, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8300, 3, 3000, 3010, 5255, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8321, 3, 3010, 3020, 5270, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8342, 3, 3020, 3030, 5285, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8363, 3, 3030, 3040, 5300, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8384, 3, 3040, 3050, 5315, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8405, 3, 3050, 3060, 5330, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8426, 3, 3080, 3090, 5360, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8447, 3, 3090, 3100, 5375, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8468, 3, 3100, 3110, 5390, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8489, 3, 3110, 3120, 5405, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8510, 3, 3120, 3130, 5420, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8531, 3, 3130, 3140, 5435, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 8552, 0, 3, 5270, 8342, 5495, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8615, 0, 3, 5285, 8363, 5540, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8678, 0, 3, 5300, 8384, 5585, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8741, 0, 3, 5315, 8405, 5630, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8804, 0, 3, 5375, 8468, 5720, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8867, 0, 3, 5390, 8489, 5765, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8930, 0, 3, 5405, 8510, 5810, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8993, 0, 3, 5420, 8531, 5855, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 9056, 0, 3, 5450, 8552, 3520, 3580,
                                                 5990, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9182, 0, 3, 5495, 8615, 3580, 3640,
                                                 6080, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9308, 0, 3, 5540, 8678, 3640, 3700,
                                                 6170, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9434, 0, 3, 5585, 8741, 3700, 3760,
                                                 6260, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9560, 0, 3, 5675, 8804, 3880, 3940,
                                                 6440, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9686, 0, 3, 5720, 8867, 3940, 4000,
                                                 6530, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9812, 0, 3, 5765, 8930, 4000, 4060,
                                                 6620, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9938, 0, 3, 5810, 8993, 4060, 4120,
                                                 6710, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10064, 0, 3, 5990, 9182, 4240, 4340,
                                                 7100, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10274, 0, 3, 6080, 9308, 4340, 4440,
                                                 7250, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10484, 0, 3, 6170, 9434, 4440, 4540,
                                                 7400, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10694, 0, 3, 6440, 9686, 4740, 4840,
                                                 7850, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10904, 0, 3, 6530, 9812, 4840, 4940,
                                                 8000, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 11114, 0, 3, 6620, 9938, 4940, 5040,
                                                 8150, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11324, 3, 5240, 5255, 8321, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11352, 3, 5255, 5270, 8342, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11380, 3, 5270, 5285, 8363, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11408, 3, 5285, 5300, 8384, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11436, 3, 5300, 5315, 8405, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11464, 3, 5345, 5360, 8447, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11492, 3, 5360, 5375, 8468, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11520, 3, 5375, 5390, 8489, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11548, 3, 5390, 5405, 8510, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11576, 3, 5405, 5420, 8531, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 11604, 0, 3, 8342, 11380, 8615, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11688, 0, 3, 8363, 11408, 8678, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11772, 0, 3, 8384, 11436, 8741, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11856, 0, 3, 8468, 11520, 8867, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11940, 0, 3, 8489, 11548, 8930, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 12024, 0, 3, 8510, 11576, 8993, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 12108, 0, 3, 8552, 11604, 5900, 5990,
                                                 9182, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12276, 0, 3, 8615, 11688, 5990, 6080,
                                                 9308, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12444, 0, 3, 8678, 11772, 6080, 6170,
                                                 9434, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12612, 0, 3, 8804, 11856, 6350, 6440,
                                                 9686, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12780, 0, 3, 8867, 11940, 6440, 6530,
                                                 9812, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12948, 0, 3, 8930, 12024, 6530, 6620,
                                                 9938, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13116, 0, 3, 9056, 12108, 6800, 6950,
                                                 10064, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13396, 0, 3, 9182, 12276, 6950, 7100,
                                                 10274, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13676, 0, 3, 9308, 12444, 7100, 7250,
                                                 10484, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13956, 0, 3, 9560, 12612, 7550, 7700,
                                                 10694, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 14236, 0, 3, 9686, 12780, 7700, 7850,
                                                 10904, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 14516, 0, 3, 9812, 12948, 7850, 8000,
                                                 11114, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14796, 3, 8300, 8321, 11352, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14832, 3, 8321, 8342, 11380, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14868, 3, 8342, 8363, 11408, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14904, 3, 8363, 8384, 11436, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14940, 3, 8426, 8447, 11492, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14976, 3, 8447, 8468, 11520, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 15012, 3, 8468, 8489, 11548, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 15048, 3, 8489, 8510, 11576, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 15084, 0, 3, 11352, 14832, 11604, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 15192, 0, 3, 11380, 14868, 11688, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 15300, 0, 3, 11408, 14904, 11772, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 15408, 0, 3, 11492, 14976, 11856, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 15516, 0, 3, 11520, 15012, 11940, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 15624, 0, 3, 11548, 15048, 12024, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 15732, 0, 3, 11604, 15192, 9056, 9182,
                                                 12276, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 15948, 0, 3, 11688, 15300, 9182, 9308,
                                                 12444, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 16164, 0, 3, 11856, 15516, 9560, 9686,
                                                 12780, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 16380, 0, 3, 11940, 15624, 9686, 9812,
                                                 12948, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 16596, 0, 3, 12276, 15948, 10064, 10274,
                                                 13676, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 16956, 0, 3, 12780, 16380, 10694, 10904,
                                                 14516, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 17316, 3, 11324, 11352, 14832, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 17361, 3, 11380, 11408, 14904, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 17406, 3, 11464, 11492, 14976, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 17451, 3, 11520, 11548, 15048, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 17496, 0, 3, 14796, 17316, 15084, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 17631, 0, 3, 14868, 17361, 15300, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 17766, 0, 3, 14940, 17406, 15408, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 17901, 0, 3, 15012, 17451, 15624, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 18036, 0, 3, 15192, 17631, 12108, 12276,
                                                 15948, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 18306, 0, 3, 15516, 17901, 12612, 12780,
                                                 16380, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 18576, 0, 3, 15732, 18036, 13116, 13396,
                                                 16596, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 19026, 0, 3, 16164, 18306, 13956, 14236,
                                                 16956, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 19476, 17766, 19026, 1, 45, ncols, alpha);

            simdgeo::geom_d_y(buffer, 19746, 17766, 19026, 1, 45, ncols, alpha);

            simdgeo::geom_d_z(buffer, 20016, 17766, 19026, 1, 45, ncols, alpha);

            simdgeo::geom_d_x(buffer, 20286, 17496, 18576, 1, 45, ncols, alpha);

            simdgeo::geom_d_y(buffer, 20556, 17496, 18576, 1, 45, ncols, alpha);

            simdgeo::geom_d_z(buffer, 20826, 17496, 18576, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 21096, 20286, 810, ncols);

            simdfunc::contract_primitives(buffer, 21906, 19476, 810, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 22716, 21906, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 22716, 17, nmax);

    simdtrf::transform_l_inner(buffer, 22716, 22176, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 85 * nvalues, nvalues, buffer, 22716, 17, nmax);

    simdtrf::transform_l_inner(buffer, 22716, 22446, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 170 * nvalues, nvalues, buffer, 22716, 17, nmax);

    simdtrf::transform_l_inner(buffer, 22716, 21096, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 255 * nvalues, nvalues, buffer, 22716, 17, nmax);

    simdtrf::transform_l_inner(buffer, 22716, 21366, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 340 * nvalues, nvalues, buffer, 22716, 17, nmax);

    simdtrf::transform_l_inner(buffer, 22716, 21636, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 425 * nvalues, nvalues, buffer, 22716, 17, nmax);
}

}  // namespace simdt2ceri
