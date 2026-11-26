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


#include "SimdElectronRepulsionRsRecGK.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_gk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 27115, 25810, 1080, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 7, 8, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 8, 9, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 9, 10, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 10, 11, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 11, 12, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 12, 13, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 13, 14, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 14, 15, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 15, 16, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 19, 20, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 20, 21, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 21, 22, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 22, 23, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 23, 24, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 24, 25, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 25, 26, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 26, 27, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 27, 28, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 51, 54, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 60, 63, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 63, 66, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 66, 69, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 69, 72, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 72, 75, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 75, 78, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 78, 81, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 81, 84, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 358, 0, 90, 96, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 373, 0, 96, 102, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 388, 0, 102, 108, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 403, 0, 108, 114, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 418, 0, 114, 120, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 433, 0, 120, 126, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 448, 0, 126, 132, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 463, 0, 144, 150, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 478, 0, 150, 156, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 493, 0, 156, 162, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 508, 0, 162, 168, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 523, 0, 168, 174, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 538, 0, 174, 180, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 553, 0, 180, 186, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 568, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 571, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 574, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 577, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 580, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 583, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 586, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 589, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 592, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 595, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 598, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 601, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 604, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 607, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 610, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 613, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 616, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 619, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 622, 3, 8, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 631, 3, 9, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 640, 3, 10, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 649, 3, 11, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 658, 3, 12, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 667, 3, 13, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 676, 3, 14, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 685, 3, 15, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 694, 3, 16, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 703, 3, 20, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 712, 3, 21, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 721, 3, 22, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 730, 3, 23, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 739, 3, 24, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 748, 3, 25, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 757, 3, 26, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 766, 3, 27, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 775, 3, 28, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 784, 0, 3, 30, 622, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 802, 0, 3, 33, 631, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 820, 0, 3, 36, 640, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 838, 0, 3, 39, 649, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 856, 0, 3, 42, 658, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 874, 0, 3, 45, 667, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 892, 0, 3, 48, 676, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 910, 0, 3, 51, 685, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 928, 0, 3, 54, 694, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 946, 0, 3, 60, 703, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 964, 0, 3, 63, 712, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 982, 0, 3, 66, 721, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1000, 0, 3, 69, 730, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1018, 0, 3, 72, 739, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1036, 0, 3, 75, 748, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1054, 0, 3, 78, 757, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1072, 0, 3, 81, 766, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1090, 0, 3, 84, 775, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1108, 0, 3, 96, 820, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1138, 0, 3, 102, 838, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1168, 0, 3, 108, 856, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1198, 0, 3, 114, 874, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1228, 0, 3, 120, 892, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1258, 0, 3, 126, 910, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1288, 0, 3, 132, 928, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1318, 0, 3, 150, 982, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1348, 0, 3, 156, 1000, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1378, 0, 3, 162, 1018, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1408, 0, 3, 168, 1036, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1438, 0, 3, 174, 1054, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1468, 0, 3, 180, 1072, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1498, 0, 3, 186, 1090, 348, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1528, 0, 3, 198, 1108, 358, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1573, 0, 3, 208, 1138, 373, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1618, 0, 3, 218, 1168, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1663, 0, 3, 228, 1198, 403, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1708, 0, 3, 238, 1228, 418, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1753, 0, 3, 248, 1258, 433, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1798, 0, 3, 258, 1288, 448, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1843, 0, 3, 278, 1318, 463, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1888, 0, 3, 288, 1348, 478, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1933, 0, 3, 298, 1378, 493, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1978, 0, 3, 308, 1408, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2023, 0, 3, 318, 1438, 523, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2068, 0, 3, 328, 1468, 538, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2113, 0, 3, 338, 1498, 553, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2158, 3, 8, 9, 571, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2164, 3, 9, 10, 574, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2170, 3, 10, 11, 577, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2176, 3, 11, 12, 580, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2182, 3, 12, 13, 583, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2188, 3, 13, 14, 586, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2194, 3, 14, 15, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2200, 3, 15, 16, 592, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2206, 3, 20, 21, 598, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2212, 3, 21, 22, 601, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2218, 3, 22, 23, 604, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2224, 3, 23, 24, 607, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2230, 3, 24, 25, 610, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2236, 3, 25, 26, 613, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2242, 3, 26, 27, 616, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2248, 3, 27, 28, 619, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2254, 0, 3, 568, 2158, 631, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2272, 0, 3, 571, 2164, 640, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2290, 0, 3, 574, 2170, 649, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2308, 0, 3, 577, 2176, 658, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2326, 0, 3, 580, 2182, 667, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2344, 0, 3, 583, 2188, 676, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2362, 0, 3, 586, 2194, 685, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2380, 0, 3, 589, 2200, 694, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2398, 0, 3, 595, 2206, 712, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2416, 0, 3, 598, 2212, 721, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2434, 0, 3, 601, 2218, 730, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2452, 0, 3, 604, 2224, 739, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2470, 0, 3, 607, 2230, 748, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2488, 0, 3, 610, 2236, 757, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2506, 0, 3, 613, 2242, 766, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2524, 0, 3, 616, 2248, 775, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2542, 0, 3, 631, 2272, 90, 96, 820,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2578, 0, 3, 640, 2290, 96, 102, 838,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2614, 0, 3, 649, 2308, 102, 108, 856,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2650, 0, 3, 658, 2326, 108, 114, 874,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2686, 0, 3, 667, 2344, 114, 120, 892,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2722, 0, 3, 676, 2362, 120, 126, 910,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2758, 0, 3, 685, 2380, 126, 132, 928,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2794, 0, 3, 712, 2416, 144, 150, 982,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2830, 0, 3, 721, 2434, 150, 156, 1000,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2866, 0, 3, 730, 2452, 156, 162, 1018,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2902, 0, 3, 739, 2470, 162, 168, 1036,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2938, 0, 3, 748, 2488, 168, 174, 1054,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2974, 0, 3, 757, 2506, 174, 180, 1072,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3010, 0, 3, 766, 2524, 180, 186, 1090,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3046, 0, 3, 820, 2578, 198, 208, 1138,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3106, 0, 3, 838, 2614, 208, 218, 1168,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3166, 0, 3, 856, 2650, 218, 228, 1198,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3226, 0, 3, 874, 2686, 228, 238, 1228,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3286, 0, 3, 892, 2722, 238, 248, 1258,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3346, 0, 3, 910, 2758, 248, 258, 1288,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3406, 0, 3, 982, 2830, 278, 288, 1348,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3466, 0, 3, 1000, 2866, 288, 298, 1378,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3526, 0, 3, 1018, 2902, 298, 308, 1408,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3586, 0, 3, 1036, 2938, 308, 318, 1438,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3646, 0, 3, 1054, 2974, 318, 328, 1468,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3706, 0, 3, 1072, 3010, 328, 338, 1498,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3766, 0, 3, 2542, 2578, 1138, 3106, 358,
                                                 373, 1618, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3856, 0, 3, 2578, 2614, 1168, 3166, 373,
                                                 388, 1663, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3946, 0, 3, 2614, 2650, 1198, 3226, 388,
                                                 403, 1708, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4036, 0, 3, 2650, 2686, 1228, 3286, 403,
                                                 418, 1753, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4126, 0, 3, 2686, 2722, 1258, 3346, 418,
                                                 433, 1798, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4216, 0, 3, 2794, 2830, 1348, 3466, 463,
                                                 478, 1933, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4306, 0, 3, 2830, 2866, 1378, 3526, 478,
                                                 493, 1978, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4396, 0, 3, 2866, 2902, 1408, 3586, 493,
                                                 508, 2023, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4486, 0, 3, 2902, 2938, 1438, 3646, 508,
                                                 523, 2068, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4576, 0, 3, 2938, 2974, 1468, 3706, 523,
                                                 538, 2113, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4666, 3, 568, 571, 2164, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4676, 3, 571, 574, 2170, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4686, 3, 574, 577, 2176, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4696, 3, 577, 580, 2182, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4706, 3, 580, 583, 2188, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4716, 3, 583, 586, 2194, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4726, 3, 586, 589, 2200, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4736, 3, 595, 598, 2212, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4746, 3, 598, 601, 2218, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4756, 3, 601, 604, 2224, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4766, 3, 604, 607, 2230, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4776, 3, 607, 610, 2236, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4786, 3, 610, 613, 2242, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4796, 3, 613, 616, 2248, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4806, 0, 3, 2158, 4666, 2272, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4836, 0, 3, 2164, 4676, 2290, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4866, 0, 3, 2170, 4686, 2308, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4896, 0, 3, 2176, 4696, 2326, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4926, 0, 3, 2182, 4706, 2344, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4956, 0, 3, 2188, 4716, 2362, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4986, 0, 3, 2194, 4726, 2380, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5016, 0, 3, 2206, 4736, 2416, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5046, 0, 3, 2212, 4746, 2434, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5076, 0, 3, 2218, 4756, 2452, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5106, 0, 3, 2224, 4766, 2470, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5136, 0, 3, 2230, 4776, 2488, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5166, 0, 3, 2236, 4786, 2506, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5196, 0, 3, 2242, 4796, 2524, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 5226, 0, 3, 2254, 4806, 784, 802, 2542,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5286, 0, 3, 2272, 4836, 802, 820, 2578,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5346, 0, 3, 2290, 4866, 820, 838, 2614,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5406, 0, 3, 2308, 4896, 838, 856, 2650,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5466, 0, 3, 2326, 4926, 856, 874, 2686,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5526, 0, 3, 2344, 4956, 874, 892, 2722,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5586, 0, 3, 2362, 4986, 892, 910, 2758,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5646, 0, 3, 2398, 5016, 946, 964, 2794,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5706, 0, 3, 2416, 5046, 964, 982, 2830,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5766, 0, 3, 2434, 5076, 982, 1000, 2866,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5826, 0, 3, 2452, 5106, 1000, 1018,
                                                 2902, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5886, 0, 3, 2470, 5136, 1018, 1036,
                                                 2938, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5946, 0, 3, 2488, 5166, 1036, 1054,
                                                 2974, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6006, 0, 3, 2506, 5196, 1054, 1072,
                                                 3010, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6066, 0, 3, 2578, 5346, 1108, 1138,
                                                 3106, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6166, 0, 3, 2614, 5406, 1138, 1168,
                                                 3166, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6266, 0, 3, 2650, 5466, 1168, 1198,
                                                 3226, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6366, 0, 3, 2686, 5526, 1198, 1228,
                                                 3286, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6466, 0, 3, 2722, 5586, 1228, 1258,
                                                 3346, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6566, 0, 3, 2830, 5766, 1318, 1348,
                                                 3466, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6666, 0, 3, 2866, 5826, 1348, 1378,
                                                 3526, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6766, 0, 3, 2902, 5886, 1378, 1408,
                                                 3586, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6866, 0, 3, 2938, 5946, 1408, 1438,
                                                 3646, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6966, 0, 3, 2974, 6006, 1438, 1468,
                                                 3706, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7066, 0, 3, 5226, 5286, 3046, 6066,
                                                 1528, 1573, 3766, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7216, 0, 3, 5286, 5346, 3106, 6166,
                                                 1573, 1618, 3856, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7366, 0, 3, 5346, 5406, 3166, 6266,
                                                 1618, 1663, 3946, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7516, 0, 3, 5406, 5466, 3226, 6366,
                                                 1663, 1708, 4036, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7666, 0, 3, 5466, 5526, 3286, 6466,
                                                 1708, 1753, 4126, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7816, 0, 3, 5646, 5706, 3406, 6566,
                                                 1843, 1888, 4216, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7966, 0, 3, 5706, 5766, 3466, 6666,
                                                 1888, 1933, 4306, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8116, 0, 3, 5766, 5826, 3526, 6766,
                                                 1933, 1978, 4396, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8266, 0, 3, 5826, 5886, 3586, 6866,
                                                 1978, 2023, 4486, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8416, 0, 3, 5886, 5946, 3646, 6966,
                                                 2023, 2068, 4576, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8566, 3, 2158, 2164, 4676, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8581, 3, 2164, 2170, 4686, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8596, 3, 2170, 2176, 4696, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8611, 3, 2176, 2182, 4706, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8626, 3, 2182, 2188, 4716, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8641, 3, 2188, 2194, 4726, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8656, 3, 2206, 2212, 4746, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8671, 3, 2212, 2218, 4756, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8686, 3, 2218, 2224, 4766, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8701, 3, 2224, 2230, 4776, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8716, 3, 2230, 2236, 4786, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8731, 3, 2236, 2242, 4796, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 8746, 0, 3, 4666, 8566, 4836, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8791, 0, 3, 4676, 8581, 4866, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8836, 0, 3, 4686, 8596, 4896, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8881, 0, 3, 4696, 8611, 4926, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8926, 0, 3, 4706, 8626, 4956, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8971, 0, 3, 4716, 8641, 4986, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9016, 0, 3, 4736, 8656, 5046, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9061, 0, 3, 4746, 8671, 5076, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9106, 0, 3, 4756, 8686, 5106, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9151, 0, 3, 4766, 8701, 5136, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9196, 0, 3, 4776, 8716, 5166, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9241, 0, 3, 4786, 8731, 5196, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 9286, 0, 3, 4836, 8791, 2542, 2578,
                                                 5346, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9376, 0, 3, 4866, 8836, 2578, 2614,
                                                 5406, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9466, 0, 3, 4896, 8881, 2614, 2650,
                                                 5466, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9556, 0, 3, 4926, 8926, 2650, 2686,
                                                 5526, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9646, 0, 3, 4956, 8971, 2686, 2722,
                                                 5586, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9736, 0, 3, 5046, 9061, 2794, 2830,
                                                 5766, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9826, 0, 3, 5076, 9106, 2830, 2866,
                                                 5826, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 9916, 0, 3, 5106, 9151, 2866, 2902,
                                                 5886, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10006, 0, 3, 5136, 9196, 2902, 2938,
                                                 5946, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10096, 0, 3, 5166, 9241, 2938, 2974,
                                                 6006, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10186, 0, 3, 5346, 9376, 3046, 3106,
                                                 6166, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10336, 0, 3, 5406, 9466, 3106, 3166,
                                                 6266, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10486, 0, 3, 5466, 9556, 3166, 3226,
                                                 6366, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10636, 0, 3, 5526, 9646, 3226, 3286,
                                                 6466, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10786, 0, 3, 5766, 9826, 3406, 3466,
                                                 6666, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10936, 0, 3, 5826, 9916, 3466, 3526,
                                                 6766, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11086, 0, 3, 5886, 10006, 3526, 3586,
                                                 6866, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11236, 0, 3, 5946, 10096, 3586, 3646,
                                                 6966, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11386, 0, 3, 9286, 9376, 6166, 10336,
                                                 3766, 3856, 7366, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11611, 0, 3, 9376, 9466, 6266, 10486,
                                                 3856, 3946, 7516, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11836, 0, 3, 9466, 9556, 6366, 10636,
                                                 3946, 4036, 7666, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 12061, 0, 3, 9736, 9826, 6666, 10936,
                                                 4216, 4306, 8116, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 12286, 0, 3, 9826, 9916, 6766, 11086,
                                                 4306, 4396, 8266, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 12511, 0, 3, 9916, 10006, 6866, 11236,
                                                 4396, 4486, 8416, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12736, 3, 4666, 4676, 8581, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12757, 3, 4676, 4686, 8596, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12778, 3, 4686, 4696, 8611, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12799, 3, 4696, 4706, 8626, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12820, 3, 4706, 4716, 8641, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12841, 3, 4736, 4746, 8671, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12862, 3, 4746, 4756, 8686, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12883, 3, 4756, 4766, 8701, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12904, 3, 4766, 4776, 8716, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12925, 3, 4776, 4786, 8731, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 12946, 0, 3, 8566, 12736, 8791, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13009, 0, 3, 8581, 12757, 8836, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13072, 0, 3, 8596, 12778, 8881, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13135, 0, 3, 8611, 12799, 8926, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13198, 0, 3, 8626, 12820, 8971, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13261, 0, 3, 8656, 12841, 9061, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13324, 0, 3, 8671, 12862, 9106, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13387, 0, 3, 8686, 12883, 9151, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13450, 0, 3, 8701, 12904, 9196, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13513, 0, 3, 8716, 12925, 9241, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 13576, 0, 3, 8746, 12946, 5226, 5286,
                                                 9286, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13702, 0, 3, 8791, 13009, 5286, 5346,
                                                 9376, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13828, 0, 3, 8836, 13072, 5346, 5406,
                                                 9466, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13954, 0, 3, 8881, 13135, 5406, 5466,
                                                 9556, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14080, 0, 3, 8926, 13198, 5466, 5526,
                                                 9646, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14206, 0, 3, 9016, 13261, 5646, 5706,
                                                 9736, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14332, 0, 3, 9061, 13324, 5706, 5766,
                                                 9826, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14458, 0, 3, 9106, 13387, 5766, 5826,
                                                 9916, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14584, 0, 3, 9151, 13450, 5826, 5886,
                                                 10006, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14710, 0, 3, 9196, 13513, 5886, 5946,
                                                 10096, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14836, 0, 3, 9376, 13828, 6066, 6166,
                                                 10336, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15046, 0, 3, 9466, 13954, 6166, 6266,
                                                 10486, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15256, 0, 3, 9556, 14080, 6266, 6366,
                                                 10636, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15466, 0, 3, 9826, 14458, 6566, 6666,
                                                 10936, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15676, 0, 3, 9916, 14584, 6666, 6766,
                                                 11086, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15886, 0, 3, 10006, 14710, 6766, 6866,
                                                 11236, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16096, 0, 3, 13576, 13702, 10186, 14836,
                                                 7066, 7216, 11386, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16411, 0, 3, 13702, 13828, 10336, 15046,
                                                 7216, 7366, 11611, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16726, 0, 3, 13828, 13954, 10486, 15256,
                                                 7366, 7516, 11836, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 17041, 0, 3, 14206, 14332, 10786, 15466,
                                                 7816, 7966, 12061, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 17356, 0, 3, 14332, 14458, 10936, 15676,
                                                 7966, 8116, 12286, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 17671, 0, 3, 14458, 14584, 11086, 15886,
                                                 8116, 8266, 12511, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17986, 3, 8566, 8581, 12757, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18014, 3, 8581, 8596, 12778, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18042, 3, 8596, 8611, 12799, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18070, 3, 8611, 8626, 12820, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18098, 3, 8656, 8671, 12862, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18126, 3, 8671, 8686, 12883, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18154, 3, 8686, 8701, 12904, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 18182, 3, 8701, 8716, 12925, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 18210, 0, 3, 12736, 17986, 13009, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18294, 0, 3, 12757, 18014, 13072, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18378, 0, 3, 12778, 18042, 13135, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18462, 0, 3, 12799, 18070, 13198, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18546, 0, 3, 12841, 18098, 13324, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18630, 0, 3, 12862, 18126, 13387, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18714, 0, 3, 12883, 18154, 13450, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18798, 0, 3, 12904, 18182, 13513, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 18882, 0, 3, 13009, 18294, 9286, 9376,
                                                 13828, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19050, 0, 3, 13072, 18378, 9376, 9466,
                                                 13954, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19218, 0, 3, 13135, 18462, 9466, 9556,
                                                 14080, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19386, 0, 3, 13324, 18630, 9736, 9826,
                                                 14458, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19554, 0, 3, 13387, 18714, 9826, 9916,
                                                 14584, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 19722, 0, 3, 13450, 18798, 9916, 10006,
                                                 14710, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 19890, 0, 3, 13828, 19050, 10186, 10336,
                                                 15046, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20170, 0, 3, 13954, 19218, 10336, 10486,
                                                 15256, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20450, 0, 3, 14458, 19554, 10786, 10936,
                                                 15676, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20730, 0, 3, 14584, 19722, 10936, 11086,
                                                 15886, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 21010, 0, 3, 18882, 19050, 15046, 20170,
                                                 11386, 11611, 16726, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 21430, 0, 3, 19386, 19554, 15676, 20730,
                                                 12061, 12286, 17671, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 21850, 3, 12736, 12757, 18014, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 21886, 3, 12757, 12778, 18042, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 21922, 3, 12778, 12799, 18070, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 21958, 3, 12841, 12862, 18126, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 21994, 3, 12862, 12883, 18154, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 22030, 3, 12883, 12904, 18182, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 22066, 0, 3, 17986, 21850, 18294, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 22174, 0, 3, 18014, 21886, 18378, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 22282, 0, 3, 18042, 21922, 18462, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 22390, 0, 3, 18098, 21958, 18630, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 22498, 0, 3, 18126, 21994, 18714, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 22606, 0, 3, 18154, 22030, 18798, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 22714, 0, 3, 18210, 22066, 13576, 13702,
                                                 18882, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 22930, 0, 3, 18294, 22174, 13702, 13828,
                                                 19050, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 23146, 0, 3, 18378, 22282, 13828, 13954,
                                                 19218, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 23362, 0, 3, 18546, 22390, 14206, 14332,
                                                 19386, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 23578, 0, 3, 18630, 22498, 14332, 14458,
                                                 19554, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 23794, 0, 3, 18714, 22606, 14458, 14584,
                                                 19722, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 24010, 0, 3, 19050, 23146, 14836, 15046,
                                                 20170, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 24370, 0, 3, 19554, 23794, 15466, 15676,
                                                 20730, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 24730, 0, 3, 22714, 22930, 19890, 24010,
                                                 16096, 16411, 21010, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 25270, 0, 3, 23362, 23578, 20450, 24370,
                                                 17041, 17356, 21430, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 25810, 24730, 1080, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 26890, 26350, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 26890, 15, nmax);

    simdtrf::transform_k_inner(buffer, 26890, 25810, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 135 * nvalues, nvalues, buffer, 26890, 15, nmax);
}

}  // namespace simdt2ceri
