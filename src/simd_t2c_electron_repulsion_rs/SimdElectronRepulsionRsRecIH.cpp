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


#include "SimdElectronRepulsionRsRecIH.hpp"

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
compute_rs_ih_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ih_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 30586, 29102, 1176, nvalues);

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

            compute_prim_hs_electron_repulsion_0(buffer, 568, 0, 198, 208, 373, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 589, 0, 208, 218, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 610, 0, 218, 228, 403, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 631, 0, 228, 238, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 652, 0, 238, 248, 433, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 673, 0, 248, 258, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 694, 0, 278, 288, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 715, 0, 288, 298, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 736, 0, 298, 308, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 757, 0, 308, 318, 523, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 778, 0, 318, 328, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 799, 0, 328, 338, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 820, 0, 358, 373, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 848, 0, 373, 388, 610, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 876, 0, 388, 403, 631, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 904, 0, 403, 418, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 932, 0, 418, 433, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 960, 0, 463, 478, 715, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 988, 0, 478, 493, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1016, 0, 493, 508, 757, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1044, 0, 508, 523, 778, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1072, 0, 523, 538, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1100, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1103, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1106, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1109, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1112, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1115, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1118, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1121, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1124, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1127, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1130, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1133, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1136, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1139, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1142, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1145, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1148, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1151, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1154, 3, 8, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1163, 3, 9, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1172, 3, 10, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1181, 3, 11, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1190, 3, 12, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1199, 3, 13, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1208, 3, 14, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1217, 3, 15, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1226, 3, 16, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1235, 3, 20, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1244, 3, 21, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1253, 3, 22, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1262, 3, 23, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1271, 3, 24, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1280, 3, 25, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1289, 3, 26, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1298, 3, 27, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1307, 3, 28, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1316, 0, 3, 30, 1154, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1334, 0, 3, 33, 1163, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1352, 0, 3, 36, 1172, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1370, 0, 3, 39, 1181, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1388, 0, 3, 42, 1190, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1406, 0, 3, 45, 1199, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1424, 0, 3, 48, 1208, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1442, 0, 3, 51, 1217, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1460, 0, 3, 54, 1226, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1478, 0, 3, 60, 1235, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1496, 0, 3, 63, 1244, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1514, 0, 3, 66, 1253, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1532, 0, 3, 69, 1262, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1550, 0, 3, 72, 1271, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1568, 0, 3, 75, 1280, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1586, 0, 3, 78, 1289, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1604, 0, 3, 81, 1298, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1622, 0, 3, 84, 1307, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1640, 0, 3, 96, 1352, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1670, 0, 3, 102, 1370, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1700, 0, 3, 108, 1388, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1730, 0, 3, 114, 1406, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1760, 0, 3, 120, 1424, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1790, 0, 3, 126, 1442, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1820, 0, 3, 132, 1460, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1850, 0, 3, 150, 1514, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1880, 0, 3, 156, 1532, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1910, 0, 3, 162, 1550, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1940, 0, 3, 168, 1568, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1970, 0, 3, 174, 1586, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2000, 0, 3, 180, 1604, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2030, 0, 3, 186, 1622, 348, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2060, 0, 3, 198, 1640, 358, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2105, 0, 3, 208, 1670, 373, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2150, 0, 3, 218, 1700, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2195, 0, 3, 228, 1730, 403, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2240, 0, 3, 238, 1760, 418, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2285, 0, 3, 248, 1790, 433, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2330, 0, 3, 258, 1820, 448, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2375, 0, 3, 278, 1850, 463, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2420, 0, 3, 288, 1880, 478, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2465, 0, 3, 298, 1910, 493, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2510, 0, 3, 308, 1940, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2555, 0, 3, 318, 1970, 523, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2600, 0, 3, 328, 2000, 538, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2645, 0, 3, 338, 2030, 553, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2690, 0, 3, 373, 2150, 589, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2753, 0, 3, 388, 2195, 610, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2816, 0, 3, 403, 2240, 631, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2879, 0, 3, 418, 2285, 652, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2942, 0, 3, 433, 2330, 673, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3005, 0, 3, 478, 2465, 715, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3068, 0, 3, 493, 2510, 736, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3131, 0, 3, 508, 2555, 757, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3194, 0, 3, 523, 2600, 778, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3257, 0, 3, 538, 2645, 799, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3320, 0, 3, 568, 2690, 820, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3404, 0, 3, 589, 2753, 848, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3488, 0, 3, 610, 2816, 876, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3572, 0, 3, 631, 2879, 904, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3656, 0, 3, 652, 2942, 932, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3740, 0, 3, 694, 3005, 960, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3824, 0, 3, 715, 3068, 988, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3908, 0, 3, 736, 3131, 1016, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3992, 0, 3, 757, 3194, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4076, 0, 3, 778, 3257, 1072, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4160, 3, 8, 9, 1103, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4166, 3, 9, 10, 1106, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4172, 3, 10, 11, 1109, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4178, 3, 11, 12, 1112, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4184, 3, 12, 13, 1115, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4190, 3, 13, 14, 1118, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4196, 3, 14, 15, 1121, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4202, 3, 15, 16, 1124, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4208, 3, 20, 21, 1130, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4214, 3, 21, 22, 1133, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4220, 3, 22, 23, 1136, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4226, 3, 23, 24, 1139, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4232, 3, 24, 25, 1142, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4238, 3, 25, 26, 1145, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4244, 3, 26, 27, 1148, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4250, 3, 27, 28, 1151, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4256, 0, 3, 1100, 4160, 1163, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4274, 0, 3, 1103, 4166, 1172, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4292, 0, 3, 1106, 4172, 1181, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4310, 0, 3, 1109, 4178, 1190, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4328, 0, 3, 1112, 4184, 1199, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4346, 0, 3, 1115, 4190, 1208, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4364, 0, 3, 1118, 4196, 1217, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4382, 0, 3, 1121, 4202, 1226, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4400, 0, 3, 1127, 4208, 1244, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4418, 0, 3, 1130, 4214, 1253, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4436, 0, 3, 1133, 4220, 1262, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4454, 0, 3, 1136, 4226, 1271, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4472, 0, 3, 1139, 4232, 1280, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4490, 0, 3, 1142, 4238, 1289, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4508, 0, 3, 1145, 4244, 1298, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4526, 0, 3, 1148, 4250, 1307, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4544, 0, 3, 1163, 4274, 90, 96, 1352,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4580, 0, 3, 1172, 4292, 96, 102, 1370,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4616, 0, 3, 1181, 4310, 102, 108, 1388,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4652, 0, 3, 1190, 4328, 108, 114, 1406,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4688, 0, 3, 1199, 4346, 114, 120, 1424,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4724, 0, 3, 1208, 4364, 120, 126, 1442,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4760, 0, 3, 1217, 4382, 126, 132, 1460,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4796, 0, 3, 1244, 4418, 144, 150, 1514,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4832, 0, 3, 1253, 4436, 150, 156, 1532,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4868, 0, 3, 1262, 4454, 156, 162, 1550,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4904, 0, 3, 1271, 4472, 162, 168, 1568,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4940, 0, 3, 1280, 4490, 168, 174, 1586,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4976, 0, 3, 1289, 4508, 174, 180, 1604,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5012, 0, 3, 1298, 4526, 180, 186, 1622,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5048, 0, 3, 1352, 4580, 198, 208, 1670,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5108, 0, 3, 1370, 4616, 208, 218, 1700,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5168, 0, 3, 1388, 4652, 218, 228, 1730,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5228, 0, 3, 1406, 4688, 228, 238, 1760,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5288, 0, 3, 1424, 4724, 238, 248, 1790,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5348, 0, 3, 1442, 4760, 248, 258, 1820,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5408, 0, 3, 1514, 4832, 278, 288, 1880,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5468, 0, 3, 1532, 4868, 288, 298, 1910,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5528, 0, 3, 1550, 4904, 298, 308, 1940,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5588, 0, 3, 1568, 4940, 308, 318, 1970,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5648, 0, 3, 1586, 4976, 318, 328, 2000,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5708, 0, 3, 1604, 5012, 328, 338, 2030,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5768, 0, 3, 4544, 4580, 1670, 5108, 358,
                                                 373, 2150, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5858, 0, 3, 4580, 4616, 1700, 5168, 373,
                                                 388, 2195, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5948, 0, 3, 4616, 4652, 1730, 5228, 388,
                                                 403, 2240, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6038, 0, 3, 4652, 4688, 1760, 5288, 403,
                                                 418, 2285, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6128, 0, 3, 4688, 4724, 1790, 5348, 418,
                                                 433, 2330, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6218, 0, 3, 4796, 4832, 1880, 5468, 463,
                                                 478, 2465, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6308, 0, 3, 4832, 4868, 1910, 5528, 478,
                                                 493, 2510, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6398, 0, 3, 4868, 4904, 1940, 5588, 493,
                                                 508, 2555, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6488, 0, 3, 4904, 4940, 1970, 5648, 508,
                                                 523, 2600, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6578, 0, 3, 4940, 4976, 2000, 5708, 523,
                                                 538, 2645, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6668, 0, 3, 5048, 5108, 2150, 5858, 568,
                                                 589, 2753, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6794, 0, 3, 5108, 5168, 2195, 5948, 589,
                                                 610, 2816, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6920, 0, 3, 5168, 5228, 2240, 6038, 610,
                                                 631, 2879, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7046, 0, 3, 5228, 5288, 2285, 6128, 631,
                                                 652, 2942, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7172, 0, 3, 5408, 5468, 2465, 6308, 694,
                                                 715, 3068, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7298, 0, 3, 5468, 5528, 2510, 6398, 715,
                                                 736, 3131, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7424, 0, 3, 5528, 5588, 2555, 6488, 736,
                                                 757, 3194, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7550, 0, 3, 5588, 5648, 2600, 6578, 757,
                                                 778, 3257, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7676, 0, 3, 5768, 5858, 2753, 6794, 820,
                                                 848, 3488, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7844, 0, 3, 5858, 5948, 2816, 6920, 848,
                                                 876, 3572, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8012, 0, 3, 5948, 6038, 2879, 7046, 876,
                                                 904, 3656, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8180, 0, 3, 6218, 6308, 3068, 7298, 960,
                                                 988, 3908, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8348, 0, 3, 6308, 6398, 3131, 7424, 988,
                                                 1016, 3992, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8516, 0, 3, 6398, 6488, 3194, 7550,
                                                 1016, 1044, 4076, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8684, 3, 1100, 1103, 4166, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8694, 3, 1103, 1106, 4172, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8704, 3, 1106, 1109, 4178, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8714, 3, 1109, 1112, 4184, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8724, 3, 1112, 1115, 4190, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8734, 3, 1115, 1118, 4196, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8744, 3, 1118, 1121, 4202, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8754, 3, 1127, 1130, 4214, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8764, 3, 1130, 1133, 4220, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8774, 3, 1133, 1136, 4226, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8784, 3, 1136, 1139, 4232, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8794, 3, 1139, 1142, 4238, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8804, 3, 1142, 1145, 4244, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8814, 3, 1145, 1148, 4250, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8824, 0, 3, 4160, 8684, 4274, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8854, 0, 3, 4166, 8694, 4292, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8884, 0, 3, 4172, 8704, 4310, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8914, 0, 3, 4178, 8714, 4328, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8944, 0, 3, 4184, 8724, 4346, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8974, 0, 3, 4190, 8734, 4364, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9004, 0, 3, 4196, 8744, 4382, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9034, 0, 3, 4208, 8754, 4418, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9064, 0, 3, 4214, 8764, 4436, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9094, 0, 3, 4220, 8774, 4454, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9124, 0, 3, 4226, 8784, 4472, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9154, 0, 3, 4232, 8794, 4490, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9184, 0, 3, 4238, 8804, 4508, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9214, 0, 3, 4244, 8814, 4526, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 9244, 0, 3, 4256, 8824, 1316, 1334,
                                                 4544, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9304, 0, 3, 4274, 8854, 1334, 1352,
                                                 4580, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9364, 0, 3, 4292, 8884, 1352, 1370,
                                                 4616, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9424, 0, 3, 4310, 8914, 1370, 1388,
                                                 4652, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9484, 0, 3, 4328, 8944, 1388, 1406,
                                                 4688, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9544, 0, 3, 4346, 8974, 1406, 1424,
                                                 4724, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9604, 0, 3, 4364, 9004, 1424, 1442,
                                                 4760, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9664, 0, 3, 4400, 9034, 1478, 1496,
                                                 4796, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9724, 0, 3, 4418, 9064, 1496, 1514,
                                                 4832, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9784, 0, 3, 4436, 9094, 1514, 1532,
                                                 4868, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9844, 0, 3, 4454, 9124, 1532, 1550,
                                                 4904, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9904, 0, 3, 4472, 9154, 1550, 1568,
                                                 4940, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9964, 0, 3, 4490, 9184, 1568, 1586,
                                                 4976, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10024, 0, 3, 4508, 9214, 1586, 1604,
                                                 5012, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10084, 0, 3, 4580, 9364, 1640, 1670,
                                                 5108, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10184, 0, 3, 4616, 9424, 1670, 1700,
                                                 5168, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10284, 0, 3, 4652, 9484, 1700, 1730,
                                                 5228, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10384, 0, 3, 4688, 9544, 1730, 1760,
                                                 5288, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10484, 0, 3, 4724, 9604, 1760, 1790,
                                                 5348, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10584, 0, 3, 4832, 9784, 1850, 1880,
                                                 5468, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10684, 0, 3, 4868, 9844, 1880, 1910,
                                                 5528, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10784, 0, 3, 4904, 9904, 1910, 1940,
                                                 5588, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10884, 0, 3, 4940, 9964, 1940, 1970,
                                                 5648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10984, 0, 3, 4976, 10024, 1970, 2000,
                                                 5708, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11084, 0, 3, 9244, 9304, 5048, 10084,
                                                 2060, 2105, 5768, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11234, 0, 3, 9304, 9364, 5108, 10184,
                                                 2105, 2150, 5858, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11384, 0, 3, 9364, 9424, 5168, 10284,
                                                 2150, 2195, 5948, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11534, 0, 3, 9424, 9484, 5228, 10384,
                                                 2195, 2240, 6038, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11684, 0, 3, 9484, 9544, 5288, 10484,
                                                 2240, 2285, 6128, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11834, 0, 3, 9664, 9724, 5408, 10584,
                                                 2375, 2420, 6218, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11984, 0, 3, 9724, 9784, 5468, 10684,
                                                 2420, 2465, 6308, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12134, 0, 3, 9784, 9844, 5528, 10784,
                                                 2465, 2510, 6398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12284, 0, 3, 9844, 9904, 5588, 10884,
                                                 2510, 2555, 6488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12434, 0, 3, 9904, 9964, 5648, 10984,
                                                 2555, 2600, 6578, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12584, 0, 3, 10084, 10184, 5858, 11384,
                                                 2690, 2753, 6794, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12794, 0, 3, 10184, 10284, 5948, 11534,
                                                 2753, 2816, 6920, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13004, 0, 3, 10284, 10384, 6038, 11684,
                                                 2816, 2879, 7046, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13214, 0, 3, 10584, 10684, 6308, 12134,
                                                 3005, 3068, 7298, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13424, 0, 3, 10684, 10784, 6398, 12284,
                                                 3068, 3131, 7424, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13634, 0, 3, 10784, 10884, 6488, 12434,
                                                 3131, 3194, 7550, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13844, 0, 3, 11084, 11234, 6668, 12584,
                                                 3320, 3404, 7676, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14124, 0, 3, 11234, 11384, 6794, 12794,
                                                 3404, 3488, 7844, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14404, 0, 3, 11384, 11534, 6920, 13004,
                                                 3488, 3572, 8012, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14684, 0, 3, 11834, 11984, 7172, 13214,
                                                 3740, 3824, 8180, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14964, 0, 3, 11984, 12134, 7298, 13424,
                                                 3824, 3908, 8348, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15244, 0, 3, 12134, 12284, 7424, 13634,
                                                 3908, 3992, 8516, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15524, 3, 4160, 4166, 8694, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15539, 3, 4166, 4172, 8704, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15554, 3, 4172, 4178, 8714, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15569, 3, 4178, 4184, 8724, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15584, 3, 4184, 4190, 8734, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15599, 3, 4190, 4196, 8744, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15614, 3, 4208, 4214, 8764, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15629, 3, 4214, 4220, 8774, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15644, 3, 4220, 4226, 8784, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15659, 3, 4226, 4232, 8794, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15674, 3, 4232, 4238, 8804, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15689, 3, 4238, 4244, 8814, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 15704, 0, 3, 8684, 15524, 8854, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15749, 0, 3, 8694, 15539, 8884, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15794, 0, 3, 8704, 15554, 8914, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15839, 0, 3, 8714, 15569, 8944, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15884, 0, 3, 8724, 15584, 8974, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15929, 0, 3, 8734, 15599, 9004, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 15974, 0, 3, 8754, 15614, 9064, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16019, 0, 3, 8764, 15629, 9094, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16064, 0, 3, 8774, 15644, 9124, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16109, 0, 3, 8784, 15659, 9154, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16154, 0, 3, 8794, 15674, 9184, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16199, 0, 3, 8804, 15689, 9214, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 16244, 0, 3, 8854, 15749, 4544, 4580,
                                                 9364, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16334, 0, 3, 8884, 15794, 4580, 4616,
                                                 9424, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16424, 0, 3, 8914, 15839, 4616, 4652,
                                                 9484, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16514, 0, 3, 8944, 15884, 4652, 4688,
                                                 9544, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16604, 0, 3, 8974, 15929, 4688, 4724,
                                                 9604, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16694, 0, 3, 9064, 16019, 4796, 4832,
                                                 9784, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16784, 0, 3, 9094, 16064, 4832, 4868,
                                                 9844, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16874, 0, 3, 9124, 16109, 4868, 4904,
                                                 9904, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16964, 0, 3, 9154, 16154, 4904, 4940,
                                                 9964, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17054, 0, 3, 9184, 16199, 4940, 4976,
                                                 10024, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17144, 0, 3, 9364, 16334, 5048, 5108,
                                                 10184, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17294, 0, 3, 9424, 16424, 5108, 5168,
                                                 10284, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17444, 0, 3, 9484, 16514, 5168, 5228,
                                                 10384, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17594, 0, 3, 9544, 16604, 5228, 5288,
                                                 10484, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17744, 0, 3, 9784, 16784, 5408, 5468,
                                                 10684, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17894, 0, 3, 9844, 16874, 5468, 5528,
                                                 10784, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18044, 0, 3, 9904, 16964, 5528, 5588,
                                                 10884, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18194, 0, 3, 9964, 17054, 5588, 5648,
                                                 10984, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18344, 0, 3, 16244, 16334, 10184, 17294,
                                                 5768, 5858, 11384, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18569, 0, 3, 16334, 16424, 10284, 17444,
                                                 5858, 5948, 11534, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18794, 0, 3, 16424, 16514, 10384, 17594,
                                                 5948, 6038, 11684, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19019, 0, 3, 16694, 16784, 10684, 17894,
                                                 6218, 6308, 12134, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19244, 0, 3, 16784, 16874, 10784, 18044,
                                                 6308, 6398, 12284, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19469, 0, 3, 16874, 16964, 10884, 18194,
                                                 6398, 6488, 12434, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 19694, 0, 3, 17144, 17294, 11384, 18569,
                                                 6668, 6794, 12794, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20009, 0, 3, 17294, 17444, 11534, 18794,
                                                 6794, 6920, 13004, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20324, 0, 3, 17744, 17894, 12134, 19244,
                                                 7172, 7298, 13424, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20639, 0, 3, 17894, 18044, 12284, 19469,
                                                 7298, 7424, 13634, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 20954, 0, 3, 18344, 18569, 12794, 20009,
                                                 7676, 7844, 14404, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 21374, 0, 3, 19019, 19244, 13424, 20639,
                                                 8180, 8348, 15244, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21794, 3, 8684, 8694, 15539, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21815, 3, 8694, 8704, 15554, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21836, 3, 8704, 8714, 15569, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21857, 3, 8714, 8724, 15584, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21878, 3, 8724, 8734, 15599, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21899, 3, 8754, 8764, 15629, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21920, 3, 8764, 8774, 15644, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21941, 3, 8774, 8784, 15659, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21962, 3, 8784, 8794, 15674, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21983, 3, 8794, 8804, 15689, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 22004, 0, 3, 15524, 21794, 15749, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22067, 0, 3, 15539, 21815, 15794, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22130, 0, 3, 15554, 21836, 15839, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22193, 0, 3, 15569, 21857, 15884, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22256, 0, 3, 15584, 21878, 15929, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22319, 0, 3, 15614, 21899, 16019, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22382, 0, 3, 15629, 21920, 16064, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22445, 0, 3, 15644, 21941, 16109, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22508, 0, 3, 15659, 21962, 16154, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 22571, 0, 3, 15674, 21983, 16199, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 22634, 0, 3, 15704, 22004, 9244, 9304,
                                                 16244, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22760, 0, 3, 15749, 22067, 9304, 9364,
                                                 16334, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22886, 0, 3, 15794, 22130, 9364, 9424,
                                                 16424, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23012, 0, 3, 15839, 22193, 9424, 9484,
                                                 16514, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23138, 0, 3, 15884, 22256, 9484, 9544,
                                                 16604, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23264, 0, 3, 15974, 22319, 9664, 9724,
                                                 16694, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23390, 0, 3, 16019, 22382, 9724, 9784,
                                                 16784, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23516, 0, 3, 16064, 22445, 9784, 9844,
                                                 16874, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23642, 0, 3, 16109, 22508, 9844, 9904,
                                                 16964, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23768, 0, 3, 16154, 22571, 9904, 9964,
                                                 17054, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23894, 0, 3, 16334, 22886, 10084, 10184,
                                                 17294, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24104, 0, 3, 16424, 23012, 10184, 10284,
                                                 17444, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24314, 0, 3, 16514, 23138, 10284, 10384,
                                                 17594, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24524, 0, 3, 16784, 23516, 10584, 10684,
                                                 17894, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24734, 0, 3, 16874, 23642, 10684, 10784,
                                                 18044, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24944, 0, 3, 16964, 23768, 10784, 10884,
                                                 18194, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25154, 0, 3, 22634, 22760, 17144, 23894,
                                                 11084, 11234, 18344, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25469, 0, 3, 22760, 22886, 17294, 24104,
                                                 11234, 11384, 18569, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25784, 0, 3, 22886, 23012, 17444, 24314,
                                                 11384, 11534, 18794, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26099, 0, 3, 23264, 23390, 17744, 24524,
                                                 11834, 11984, 19019, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26414, 0, 3, 23390, 23516, 17894, 24734,
                                                 11984, 12134, 19244, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26729, 0, 3, 23516, 23642, 18044, 24944,
                                                 12134, 12284, 19469, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 27044, 0, 3, 23894, 24104, 18569, 25784,
                                                 12584, 12794, 20009, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 27485, 0, 3, 24524, 24734, 19244, 26729,
                                                 13214, 13424, 20639, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 27926, 0, 3, 25154, 25469, 19694, 27044,
                                                 13844, 14124, 20954, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 28514, 0, 3, 26099, 26414, 20324, 27485,
                                                 14684, 14964, 21374, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 29102, 27926, 1176, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 30278, 29690, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 30278, 11, nmax);

    simdtrf::transform_h_inner(buffer, 30278, 29102, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 143 * nvalues, nvalues, buffer, 30278, 11, nmax);
}

}  // namespace simdt2ceri
