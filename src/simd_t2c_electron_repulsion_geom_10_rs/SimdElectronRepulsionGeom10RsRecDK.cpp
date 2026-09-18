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


#include "SimdElectronRepulsionGeom10RsRecDK.hpp"

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
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
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
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_dk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 15688, 14302, 1296, nvalues);

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
                                                9, 10}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 17, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 27, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 7, 8, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 8, 9, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 9, 10, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 10, 11, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 11, 12, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 12, 13, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 13, 14, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 14, 15, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 20, 21, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 21, 22, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 22, 23, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 23, 24, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 24, 25, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 25, 26, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 28, 31, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 31, 34, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 34, 37, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 37, 40, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 40, 43, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 43, 46, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 46, 49, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 55, 58, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 58, 61, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 61, 64, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 64, 67, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 67, 70, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 70, 73, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 73, 76, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 318, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 321, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 324, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 327, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 330, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 333, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 336, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 339, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 342, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 345, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 348, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 351, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 354, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 357, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 360, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 363, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 366, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 369, 3, 27, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 372, 3, 9, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 381, 3, 10, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 390, 3, 11, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 399, 3, 12, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 408, 3, 13, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 417, 3, 14, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 426, 3, 15, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 435, 3, 20, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 444, 3, 21, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 453, 3, 22, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 462, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 471, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 480, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 489, 3, 26, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 498, 0, 3, 31, 372, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 516, 0, 3, 34, 381, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 534, 0, 3, 37, 390, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 552, 0, 3, 40, 399, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 570, 0, 3, 43, 408, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 588, 0, 3, 46, 417, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 606, 0, 3, 49, 426, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 624, 0, 3, 58, 435, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 642, 0, 3, 61, 444, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 660, 0, 3, 64, 453, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 678, 0, 3, 67, 462, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 696, 0, 3, 70, 471, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 714, 0, 3, 73, 480, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 732, 0, 3, 76, 489, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 750, 0, 3, 82, 498, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 780, 0, 3, 88, 516, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 810, 0, 3, 94, 534, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 840, 0, 3, 100, 552, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 870, 0, 3, 106, 570, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 900, 0, 3, 112, 588, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 930, 0, 3, 118, 606, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 960, 0, 3, 130, 624, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 990, 0, 3, 136, 642, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1020, 0, 3, 142, 660, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1050, 0, 3, 148, 678, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1080, 0, 3, 154, 696, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1110, 0, 3, 160, 714, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1140, 0, 3, 166, 732, 308, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1170, 3, 7, 8, 321, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1176, 3, 8, 9, 324, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1182, 3, 9, 10, 327, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1188, 3, 10, 11, 330, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1194, 3, 11, 12, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1200, 3, 12, 13, 336, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1206, 3, 13, 14, 339, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1212, 3, 14, 15, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1218, 3, 18, 19, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1224, 3, 19, 20, 351, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1230, 3, 20, 21, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1236, 3, 21, 22, 357, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1242, 3, 22, 23, 360, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1248, 3, 23, 24, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1254, 3, 24, 25, 366, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1260, 3, 25, 26, 369, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1266, 0, 3, 324, 1182, 381, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1284, 0, 3, 327, 1188, 390, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1302, 0, 3, 330, 1194, 399, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1320, 0, 3, 333, 1200, 408, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1338, 0, 3, 336, 1206, 417, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1356, 0, 3, 339, 1212, 426, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1374, 0, 3, 351, 1230, 444, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1392, 0, 3, 354, 1236, 453, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1410, 0, 3, 357, 1242, 462, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1428, 0, 3, 360, 1248, 471, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1446, 0, 3, 363, 1254, 480, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1464, 0, 3, 366, 1260, 489, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1482, 0, 3, 372, 1266, 82, 88, 516,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1518, 0, 3, 381, 1284, 88, 94, 534,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1554, 0, 3, 390, 1302, 94, 100, 552,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1590, 0, 3, 399, 1320, 100, 106, 570,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1626, 0, 3, 408, 1338, 106, 112, 588,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1662, 0, 3, 417, 1356, 112, 118, 606,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1698, 0, 3, 435, 1374, 130, 136, 642,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1734, 0, 3, 444, 1392, 136, 142, 660,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1770, 0, 3, 453, 1410, 142, 148, 678,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1806, 0, 3, 462, 1428, 148, 154, 696,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1842, 0, 3, 471, 1446, 154, 160, 714,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1878, 0, 3, 480, 1464, 160, 166, 732,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1914, 0, 3, 516, 1518, 178, 188, 810,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1974, 0, 3, 534, 1554, 188, 198, 840,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2034, 0, 3, 552, 1590, 198, 208, 870,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2094, 0, 3, 570, 1626, 208, 218, 900,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2154, 0, 3, 588, 1662, 218, 228, 930,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2214, 0, 3, 642, 1734, 248, 258, 1020,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2274, 0, 3, 660, 1770, 258, 268, 1050,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2334, 0, 3, 678, 1806, 268, 278, 1080,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2394, 0, 3, 696, 1842, 278, 288, 1110,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2454, 0, 3, 714, 1878, 288, 298, 1140,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2514, 3, 318, 321, 1176, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2524, 3, 321, 324, 1182, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2534, 3, 324, 327, 1188, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2544, 3, 327, 330, 1194, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2554, 3, 330, 333, 1200, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2564, 3, 333, 336, 1206, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2574, 3, 336, 339, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2584, 3, 345, 348, 1224, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2594, 3, 348, 351, 1230, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2604, 3, 351, 354, 1236, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2614, 3, 354, 357, 1242, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2624, 3, 357, 360, 1248, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2634, 3, 360, 363, 1254, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2644, 3, 363, 366, 1260, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2654, 0, 3, 1182, 2534, 1284, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2684, 0, 3, 1188, 2544, 1302, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2714, 0, 3, 1194, 2554, 1320, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2744, 0, 3, 1200, 2564, 1338, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2774, 0, 3, 1206, 2574, 1356, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2804, 0, 3, 1230, 2604, 1392, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2834, 0, 3, 1236, 2614, 1410, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2864, 0, 3, 1242, 2624, 1428, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2894, 0, 3, 1248, 2634, 1446, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2924, 0, 3, 1254, 2644, 1464, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2954, 0, 3, 1266, 2654, 498, 516, 1518,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3014, 0, 3, 1284, 2684, 516, 534, 1554,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3074, 0, 3, 1302, 2714, 534, 552, 1590,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3134, 0, 3, 1320, 2744, 552, 570, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3194, 0, 3, 1338, 2774, 570, 588, 1662,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3254, 0, 3, 1374, 2804, 624, 642, 1734,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3314, 0, 3, 1392, 2834, 642, 660, 1770,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3374, 0, 3, 1410, 2864, 660, 678, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3434, 0, 3, 1428, 2894, 678, 696, 1842,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3494, 0, 3, 1446, 2924, 696, 714, 1878,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3554, 0, 3, 1482, 2954, 750, 780, 1914,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3654, 0, 3, 1518, 3014, 780, 810, 1974,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3754, 0, 3, 1554, 3074, 810, 840, 2034,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3854, 0, 3, 1590, 3134, 840, 870, 2094,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3954, 0, 3, 1626, 3194, 870, 900, 2154,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4054, 0, 3, 1698, 3254, 960, 990, 2214,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4154, 0, 3, 1734, 3314, 990, 1020, 2274,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4254, 0, 3, 1770, 3374, 1020, 1050,
                                                 2334, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4354, 0, 3, 1806, 3434, 1050, 1080,
                                                 2394, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4454, 0, 3, 1842, 3494, 1080, 1110,
                                                 2454, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4554, 3, 1170, 1176, 2524, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4569, 3, 1176, 1182, 2534, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4584, 3, 1182, 1188, 2544, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4599, 3, 1188, 1194, 2554, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4614, 3, 1194, 1200, 2564, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4629, 3, 1200, 1206, 2574, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4644, 3, 1218, 1224, 2594, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4659, 3, 1224, 1230, 2604, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4674, 3, 1230, 1236, 2614, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4689, 3, 1236, 1242, 2624, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4704, 3, 1242, 1248, 2634, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4719, 3, 1248, 1254, 2644, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4734, 0, 3, 2534, 4584, 2684, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4779, 0, 3, 2544, 4599, 2714, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4824, 0, 3, 2554, 4614, 2744, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4869, 0, 3, 2564, 4629, 2774, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4914, 0, 3, 2604, 4674, 2834, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4959, 0, 3, 2614, 4689, 2864, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5004, 0, 3, 2624, 4704, 2894, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5049, 0, 3, 2634, 4719, 2924, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 5094, 0, 3, 2654, 4734, 1482, 1518,
                                                 3014, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5184, 0, 3, 2684, 4779, 1518, 1554,
                                                 3074, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5274, 0, 3, 2714, 4824, 1554, 1590,
                                                 3134, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5364, 0, 3, 2744, 4869, 1590, 1626,
                                                 3194, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5454, 0, 3, 2804, 4914, 1698, 1734,
                                                 3314, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5544, 0, 3, 2834, 4959, 1734, 1770,
                                                 3374, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5634, 0, 3, 2864, 5004, 1770, 1806,
                                                 3434, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5724, 0, 3, 2894, 5049, 1806, 1842,
                                                 3494, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5814, 0, 3, 3014, 5184, 1914, 1974,
                                                 3754, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5964, 0, 3, 3074, 5274, 1974, 2034,
                                                 3854, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6114, 0, 3, 3134, 5364, 2034, 2094,
                                                 3954, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6264, 0, 3, 3314, 5544, 2214, 2274,
                                                 4254, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6414, 0, 3, 3374, 5634, 2274, 2334,
                                                 4354, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6564, 0, 3, 3434, 5724, 2334, 2394,
                                                 4454, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6714, 3, 2514, 2524, 4569, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6735, 3, 2524, 2534, 4584, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6756, 3, 2534, 2544, 4599, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6777, 3, 2544, 2554, 4614, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6798, 3, 2554, 2564, 4629, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6819, 3, 2584, 2594, 4659, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6840, 3, 2594, 2604, 4674, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6861, 3, 2604, 2614, 4689, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6882, 3, 2614, 2624, 4704, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6903, 3, 2624, 2634, 4719, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 6924, 0, 3, 4584, 6756, 4779, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6987, 0, 3, 4599, 6777, 4824, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7050, 0, 3, 4614, 6798, 4869, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7113, 0, 3, 4674, 6861, 4959, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7176, 0, 3, 4689, 6882, 5004, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7239, 0, 3, 4704, 6903, 5049, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 7302, 0, 3, 4734, 6924, 2954, 3014,
                                                 5184, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7428, 0, 3, 4779, 6987, 3014, 3074,
                                                 5274, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7554, 0, 3, 4824, 7050, 3074, 3134,
                                                 5364, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7680, 0, 3, 4914, 7113, 3254, 3314,
                                                 5544, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7806, 0, 3, 4959, 7176, 3314, 3374,
                                                 5634, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7932, 0, 3, 5004, 7239, 3374, 3434,
                                                 5724, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8058, 0, 3, 5094, 7302, 3554, 3654,
                                                 5814, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8268, 0, 3, 5184, 7428, 3654, 3754,
                                                 5964, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8478, 0, 3, 5274, 7554, 3754, 3854,
                                                 6114, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8688, 0, 3, 5454, 7680, 4054, 4154,
                                                 6264, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8898, 0, 3, 5544, 7806, 4154, 4254,
                                                 6414, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9108, 0, 3, 5634, 7932, 4254, 4354,
                                                 6564, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9318, 3, 4554, 4569, 6735, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9346, 3, 4569, 4584, 6756, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9374, 3, 4584, 4599, 6777, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9402, 3, 4599, 4614, 6798, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9430, 3, 4644, 4659, 6840, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9458, 3, 4659, 4674, 6861, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9486, 3, 4674, 4689, 6882, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9514, 3, 4689, 4704, 6903, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 9542, 0, 3, 6735, 9346, 6924, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9626, 0, 3, 6756, 9374, 6987, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9710, 0, 3, 6777, 9402, 7050, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9794, 0, 3, 6840, 9458, 7113, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9878, 0, 3, 6861, 9486, 7176, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9962, 0, 3, 6882, 9514, 7239, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 10046, 0, 3, 6924, 9626, 5094, 5184,
                                                 7428, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 10214, 0, 3, 6987, 9710, 5184, 5274,
                                                 7554, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 10382, 0, 3, 7113, 9878, 5454, 5544,
                                                 7806, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 10550, 0, 3, 7176, 9962, 5544, 5634,
                                                 7932, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 10718, 0, 3, 7428, 10214, 5814, 5964,
                                                 8478, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 10998, 0, 3, 7806, 10550, 6264, 6414,
                                                 9108, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 11278, 3, 6714, 6735, 9346, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 11314, 3, 6756, 6777, 9402, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 11350, 3, 6819, 6840, 9458, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 11386, 3, 6861, 6882, 9514, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 11422, 0, 3, 9318, 11278, 9542, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 11530, 0, 3, 9374, 11314, 9710, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 11638, 0, 3, 9430, 11350, 9794, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 11746, 0, 3, 9486, 11386, 9962, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 11854, 0, 3, 9626, 11530, 7302, 7428,
                                                 10214, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 12070, 0, 3, 9878, 11746, 7680, 7806,
                                                 10550, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 12286, 0, 3, 10046, 11854, 8058, 8268,
                                                 10718, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 12646, 0, 3, 10382, 12070, 8688, 8898,
                                                 10998, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 13006, 11638, 12646, 1, 36, ncols, alpha);

            simdgeo::geom_d_y(buffer, 13222, 11638, 12646, 1, 36, ncols, alpha);

            simdgeo::geom_d_z(buffer, 13438, 11638, 12646, 1, 36, ncols, alpha);

            simdgeo::geom_d_x(buffer, 13654, 11422, 12286, 1, 36, ncols, alpha);

            simdgeo::geom_d_y(buffer, 13870, 11422, 12286, 1, 36, ncols, alpha);

            simdgeo::geom_d_z(buffer, 14086, 11422, 12286, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 14302, 13654, 648, ncols);

            simdfunc::contract_primitives(buffer, 14950, 13006, 648, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 15598, 14950, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 15598, 15, nmax);

    simdtrf::transform_k_inner(buffer, 15598, 15166, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 75 * nvalues, nvalues, buffer, 15598, 15, nmax);

    simdtrf::transform_k_inner(buffer, 15598, 15382, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 150 * nvalues, nvalues, buffer, 15598, 15, nmax);

    simdtrf::transform_k_inner(buffer, 15598, 14302, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 225 * nvalues, nvalues, buffer, 15598, 15, nmax);

    simdtrf::transform_k_inner(buffer, 15598, 14518, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 300 * nvalues, nvalues, buffer, 15598, 15, nmax);

    simdtrf::transform_k_inner(buffer, 15598, 14734, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 375 * nvalues, nvalues, buffer, 15598, 15, nmax);
}

}  // namespace simdt2ceri
