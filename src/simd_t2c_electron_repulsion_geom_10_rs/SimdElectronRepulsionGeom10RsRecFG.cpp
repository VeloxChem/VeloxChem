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


#include "SimdElectronRepulsionGeom10RsRecFG.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryF1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_fg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_fg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 7656, 6666, 900, nvalues);

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

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 47, 50, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 50, 53, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 53, 56, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 56, 59, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 59, 62, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 252, 0, 68, 74, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 267, 0, 74, 80, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 282, 0, 80, 86, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 297, 0, 86, 92, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 312, 0, 92, 98, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 327, 0, 110, 116, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 342, 0, 116, 122, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 357, 0, 122, 128, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 372, 0, 128, 134, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 387, 0, 134, 140, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 402, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 405, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 408, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 411, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 414, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 417, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 420, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 423, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 426, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 429, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 432, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 435, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 438, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 447, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 456, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 465, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 474, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 483, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 492, 3, 19, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 501, 3, 20, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 510, 3, 21, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 519, 3, 22, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 528, 3, 23, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 537, 3, 24, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 546, 0, 3, 29, 447, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 564, 0, 3, 32, 456, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 582, 0, 3, 35, 465, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 600, 0, 3, 38, 474, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 618, 0, 3, 41, 483, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 636, 0, 3, 50, 501, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 654, 0, 3, 53, 510, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 672, 0, 3, 56, 519, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 690, 0, 3, 59, 528, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 708, 0, 3, 62, 537, 146, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 726, 0, 3, 80, 564, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 756, 0, 3, 86, 582, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 786, 0, 3, 92, 600, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 816, 0, 3, 98, 618, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 846, 0, 3, 122, 654, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 876, 0, 3, 128, 672, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 906, 0, 3, 134, 690, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 936, 0, 3, 140, 708, 242, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 966, 0, 3, 162, 756, 282, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1011, 0, 3, 172, 786, 297, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1056, 0, 3, 182, 816, 312, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1101, 0, 3, 212, 876, 357, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1146, 0, 3, 222, 906, 372, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1191, 0, 3, 232, 936, 387, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1236, 3, 9, 10, 405, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1242, 3, 10, 11, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1248, 3, 11, 12, 411, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1254, 3, 12, 13, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1260, 3, 13, 14, 417, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1266, 3, 19, 20, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1272, 3, 20, 21, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1278, 3, 21, 22, 429, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1284, 3, 22, 23, 432, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1290, 3, 23, 24, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1296, 0, 3, 402, 1236, 447, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1314, 0, 3, 405, 1242, 456, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1332, 0, 3, 408, 1248, 465, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1350, 0, 3, 411, 1254, 474, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1368, 0, 3, 414, 1260, 483, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1386, 0, 3, 420, 1266, 501, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1404, 0, 3, 423, 1272, 510, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1422, 0, 3, 426, 1278, 519, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1440, 0, 3, 429, 1284, 528, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1458, 0, 3, 432, 1290, 537, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1476, 0, 3, 438, 1296, 68, 74, 546,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1512, 0, 3, 447, 1314, 74, 80, 564,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1548, 0, 3, 456, 1332, 80, 86, 582,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1584, 0, 3, 465, 1350, 86, 92, 600,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1620, 0, 3, 474, 1368, 92, 98, 618,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1656, 0, 3, 492, 1386, 110, 116, 636,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1692, 0, 3, 501, 1404, 116, 122, 654,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1728, 0, 3, 510, 1422, 122, 128, 672,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1764, 0, 3, 519, 1440, 128, 134, 690,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1800, 0, 3, 528, 1458, 134, 140, 708,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1836, 0, 3, 564, 1548, 152, 162, 756,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1896, 0, 3, 582, 1584, 162, 172, 786,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1956, 0, 3, 600, 1620, 172, 182, 816,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2016, 0, 3, 654, 1728, 202, 212, 876,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2076, 0, 3, 672, 1764, 212, 222, 906,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2136, 0, 3, 690, 1800, 222, 232, 936,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2196, 0, 3, 1476, 1512, 726, 1836, 252,
                                                 267, 966, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2286, 0, 3, 1512, 1548, 756, 1896, 267,
                                                 282, 1011, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2376, 0, 3, 1548, 1584, 786, 1956, 282,
                                                 297, 1056, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2466, 0, 3, 1656, 1692, 846, 2016, 327,
                                                 342, 1101, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2556, 0, 3, 1692, 1728, 876, 2076, 342,
                                                 357, 1146, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2646, 0, 3, 1728, 1764, 906, 2136, 357,
                                                 372, 1191, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2736, 3, 402, 405, 1242, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2746, 3, 405, 408, 1248, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2756, 3, 408, 411, 1254, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2766, 3, 411, 414, 1260, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2776, 3, 420, 423, 1272, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2786, 3, 423, 426, 1278, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2796, 3, 426, 429, 1284, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2806, 3, 429, 432, 1290, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2816, 0, 3, 1236, 2736, 1314, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2846, 0, 3, 1242, 2746, 1332, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2876, 0, 3, 1248, 2756, 1350, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2906, 0, 3, 1254, 2766, 1368, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2936, 0, 3, 1266, 2776, 1404, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2966, 0, 3, 1272, 2786, 1422, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2996, 0, 3, 1278, 2796, 1440, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3026, 0, 3, 1284, 2806, 1458, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3056, 0, 3, 1314, 2846, 546, 564, 1548,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3116, 0, 3, 1332, 2876, 564, 582, 1584,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3176, 0, 3, 1350, 2906, 582, 600, 1620,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3236, 0, 3, 1404, 2966, 636, 654, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3296, 0, 3, 1422, 2996, 654, 672, 1764,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3356, 0, 3, 1440, 3026, 672, 690, 1800,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3416, 0, 3, 1548, 3116, 726, 756, 1896,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3516, 0, 3, 1584, 3176, 756, 786, 1956,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3616, 0, 3, 1728, 3296, 846, 876, 2076,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3716, 0, 3, 1764, 3356, 876, 906, 2136,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3816, 0, 3, 3056, 3116, 1896, 3516, 966,
                                                 1011, 2376, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3966, 0, 3, 3236, 3296, 2076, 3716,
                                                 1101, 1146, 2646, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4116, 3, 1236, 1242, 2746, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4131, 3, 1242, 1248, 2756, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4146, 3, 1248, 1254, 2766, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4161, 3, 1266, 1272, 2786, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4176, 3, 1272, 1278, 2796, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4191, 3, 1278, 1284, 2806, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4206, 0, 3, 2736, 4116, 2846, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4251, 0, 3, 2746, 4131, 2876, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4296, 0, 3, 2756, 4146, 2906, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4341, 0, 3, 2776, 4161, 2966, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4386, 0, 3, 2786, 4176, 2996, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4431, 0, 3, 2796, 4191, 3026, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 4476, 0, 3, 2816, 4206, 1476, 1512,
                                                 3056, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4566, 0, 3, 2846, 4251, 1512, 1548,
                                                 3116, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4656, 0, 3, 2876, 4296, 1548, 1584,
                                                 3176, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4746, 0, 3, 2936, 4341, 1656, 1692,
                                                 3236, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4836, 0, 3, 2966, 4386, 1692, 1728,
                                                 3296, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4926, 0, 3, 2996, 4431, 1728, 1764,
                                                 3356, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5016, 0, 3, 3116, 4656, 1836, 1896,
                                                 3516, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5166, 0, 3, 3296, 4926, 2016, 2076,
                                                 3716, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5316, 0, 3, 4476, 4566, 3416, 5016,
                                                 2196, 2286, 3816, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5541, 0, 3, 4746, 4836, 3616, 5166,
                                                 2466, 2556, 3966, ncols, alpha, beta, p);

            simdgeo::geom_f_x(buffer, 5766, 4746, 5541, 1, 15, ncols, alpha);

            simdgeo::geom_f_y(buffer, 5916, 4746, 5541, 1, 15, ncols, alpha);

            simdgeo::geom_f_z(buffer, 6066, 4746, 5541, 1, 15, ncols, alpha);

            simdgeo::geom_f_x(buffer, 6216, 4476, 5316, 1, 15, ncols, alpha);

            simdgeo::geom_f_y(buffer, 6366, 4476, 5316, 1, 15, ncols, alpha);

            simdgeo::geom_f_z(buffer, 6516, 4476, 5316, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 6666, 6216, 450, ncols);

            simdfunc::contract_primitives(buffer, 7116, 5766, 450, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 7566, 7116, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 7566, 9, nmax);

    simdtrf::transform_g_inner(buffer, 7566, 7266, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 63 * nvalues, nvalues, buffer, 7566, 9, nmax);

    simdtrf::transform_g_inner(buffer, 7566, 7416, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 126 * nvalues, nvalues, buffer, 7566, 9, nmax);

    simdtrf::transform_g_inner(buffer, 7566, 6666, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 189 * nvalues, nvalues, buffer, 7566, 9, nmax);

    simdtrf::transform_g_inner(buffer, 7566, 6816, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 252 * nvalues, nvalues, buffer, 7566, 9, nmax);

    simdtrf::transform_g_inner(buffer, 7566, 6966, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 315 * nvalues, nvalues, buffer, 7566, 9, nmax);
}

}  // namespace simdt2ceri
