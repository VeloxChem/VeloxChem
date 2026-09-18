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


#include "SimdElectronRepulsionRsRecFL.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_fl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 19494, 18424, 900, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 384, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 387, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 390, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 393, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 396, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 399, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 402, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 405, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 408, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 411, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 414, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 417, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 420, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 423, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 426, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 429, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 432, 3, 9, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 441, 3, 10, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 450, 3, 11, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 459, 3, 12, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 468, 3, 13, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 477, 3, 14, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 486, 3, 15, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 495, 3, 16, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 504, 3, 21, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 513, 3, 22, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 522, 3, 23, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 531, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 540, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 549, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 558, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 567, 3, 28, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 576, 0, 3, 36, 432, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 594, 0, 3, 39, 441, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 612, 0, 3, 42, 450, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 630, 0, 3, 45, 459, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 648, 0, 3, 48, 468, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 666, 0, 3, 51, 477, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 684, 0, 3, 54, 486, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 702, 0, 3, 57, 495, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 720, 0, 3, 69, 504, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 738, 0, 3, 72, 513, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 756, 0, 3, 75, 522, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 774, 0, 3, 78, 531, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 792, 0, 3, 81, 540, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 810, 0, 3, 84, 549, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 828, 0, 3, 87, 558, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 846, 0, 3, 90, 567, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 864, 0, 3, 102, 594, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 894, 0, 3, 108, 612, 234, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 924, 0, 3, 114, 630, 244, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 954, 0, 3, 120, 648, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 984, 0, 3, 126, 666, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1014, 0, 3, 132, 684, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1044, 0, 3, 138, 702, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1074, 0, 3, 156, 738, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1104, 0, 3, 162, 756, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1134, 0, 3, 168, 774, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1164, 0, 3, 174, 792, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1194, 0, 3, 180, 810, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1224, 0, 3, 186, 828, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1254, 0, 3, 192, 846, 374, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1284, 3, 9, 10, 387, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1290, 3, 10, 11, 390, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1296, 3, 11, 12, 393, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1302, 3, 12, 13, 396, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1308, 3, 13, 14, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1314, 3, 14, 15, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1320, 3, 15, 16, 405, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1326, 3, 21, 22, 411, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1332, 3, 22, 23, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1338, 3, 23, 24, 417, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1344, 3, 24, 25, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1350, 3, 25, 26, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1356, 3, 26, 27, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1362, 3, 27, 28, 429, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1368, 0, 3, 384, 1284, 441, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1386, 0, 3, 387, 1290, 450, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1404, 0, 3, 390, 1296, 459, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1422, 0, 3, 393, 1302, 468, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1440, 0, 3, 396, 1308, 477, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1458, 0, 3, 399, 1314, 486, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1476, 0, 3, 402, 1320, 495, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1494, 0, 3, 408, 1326, 513, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1512, 0, 3, 411, 1332, 522, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1530, 0, 3, 414, 1338, 531, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1548, 0, 3, 417, 1344, 540, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1566, 0, 3, 420, 1350, 549, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1584, 0, 3, 423, 1356, 558, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1602, 0, 3, 426, 1362, 567, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1620, 0, 3, 432, 1368, 96, 102, 594,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1656, 0, 3, 441, 1386, 102, 108, 612,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1692, 0, 3, 450, 1404, 108, 114, 630,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1728, 0, 3, 459, 1422, 114, 120, 648,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1764, 0, 3, 468, 1440, 120, 126, 666,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1800, 0, 3, 477, 1458, 126, 132, 684,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1836, 0, 3, 486, 1476, 132, 138, 702,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1872, 0, 3, 504, 1494, 150, 156, 738,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1908, 0, 3, 513, 1512, 156, 162, 756,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1944, 0, 3, 522, 1530, 162, 168, 774,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1980, 0, 3, 531, 1548, 168, 174, 792,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2016, 0, 3, 540, 1566, 174, 180, 810,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2052, 0, 3, 549, 1584, 180, 186, 828,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2088, 0, 3, 558, 1602, 186, 192, 846,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2124, 0, 3, 576, 1620, 204, 214, 864,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2184, 0, 3, 594, 1656, 214, 224, 894,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2244, 0, 3, 612, 1692, 224, 234, 924,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2304, 0, 3, 630, 1728, 234, 244, 954,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2364, 0, 3, 648, 1764, 244, 254, 984,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2424, 0, 3, 666, 1800, 254, 264, 1014,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2484, 0, 3, 684, 1836, 264, 274, 1044,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2544, 0, 3, 720, 1872, 294, 304, 1074,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2604, 0, 3, 738, 1908, 304, 314, 1104,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2664, 0, 3, 756, 1944, 314, 324, 1134,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2724, 0, 3, 774, 1980, 324, 334, 1164,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2784, 0, 3, 792, 2016, 334, 344, 1194,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2844, 0, 3, 810, 2052, 344, 354, 1224,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2904, 0, 3, 828, 2088, 354, 364, 1254,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2964, 3, 384, 387, 1290, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2974, 3, 387, 390, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2984, 3, 390, 393, 1302, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2994, 3, 393, 396, 1308, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3004, 3, 396, 399, 1314, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3014, 3, 399, 402, 1320, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3024, 3, 408, 411, 1332, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3034, 3, 411, 414, 1338, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3044, 3, 414, 417, 1344, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3054, 3, 417, 420, 1350, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3064, 3, 420, 423, 1356, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3074, 3, 423, 426, 1362, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3084, 0, 3, 1284, 2964, 1386, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3114, 0, 3, 1290, 2974, 1404, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3144, 0, 3, 1296, 2984, 1422, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3174, 0, 3, 1302, 2994, 1440, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3204, 0, 3, 1308, 3004, 1458, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3234, 0, 3, 1314, 3014, 1476, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3264, 0, 3, 1326, 3024, 1512, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3294, 0, 3, 1332, 3034, 1530, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3324, 0, 3, 1338, 3044, 1548, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3354, 0, 3, 1344, 3054, 1566, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3384, 0, 3, 1350, 3064, 1584, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3414, 0, 3, 1356, 3074, 1602, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3444, 0, 3, 1368, 3084, 576, 594, 1656,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3504, 0, 3, 1386, 3114, 594, 612, 1692,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3564, 0, 3, 1404, 3144, 612, 630, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3624, 0, 3, 1422, 3174, 630, 648, 1764,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3684, 0, 3, 1440, 3204, 648, 666, 1800,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3744, 0, 3, 1458, 3234, 666, 684, 1836,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3804, 0, 3, 1494, 3264, 720, 738, 1908,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3864, 0, 3, 1512, 3294, 738, 756, 1944,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3924, 0, 3, 1530, 3324, 756, 774, 1980,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3984, 0, 3, 1548, 3354, 774, 792, 2016,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4044, 0, 3, 1566, 3384, 792, 810, 2052,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4104, 0, 3, 1584, 3414, 810, 828, 2088,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4164, 0, 3, 1656, 3504, 864, 894, 2244,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4264, 0, 3, 1692, 3564, 894, 924, 2304,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4364, 0, 3, 1728, 3624, 924, 954, 2364,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4464, 0, 3, 1764, 3684, 954, 984, 2424,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4564, 0, 3, 1800, 3744, 984, 1014, 2484,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4664, 0, 3, 1908, 3864, 1074, 1104,
                                                 2664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4764, 0, 3, 1944, 3924, 1104, 1134,
                                                 2724, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4864, 0, 3, 1980, 3984, 1134, 1164,
                                                 2784, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4964, 0, 3, 2016, 4044, 1164, 1194,
                                                 2844, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5064, 0, 3, 2052, 4104, 1194, 1224,
                                                 2904, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5164, 3, 1284, 1290, 2974, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5179, 3, 1290, 1296, 2984, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5194, 3, 1296, 1302, 2994, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5209, 3, 1302, 1308, 3004, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5224, 3, 1308, 1314, 3014, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5239, 3, 1326, 1332, 3034, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5254, 3, 1332, 1338, 3044, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5269, 3, 1338, 1344, 3054, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5284, 3, 1344, 1350, 3064, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5299, 3, 1350, 1356, 3074, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 5314, 0, 3, 2964, 5164, 3114, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5359, 0, 3, 2974, 5179, 3144, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5404, 0, 3, 2984, 5194, 3174, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5449, 0, 3, 2994, 5209, 3204, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5494, 0, 3, 3004, 5224, 3234, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5539, 0, 3, 3024, 5239, 3294, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5584, 0, 3, 3034, 5254, 3324, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5629, 0, 3, 3044, 5269, 3354, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5674, 0, 3, 3054, 5284, 3384, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5719, 0, 3, 3064, 5299, 3414, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 5764, 0, 3, 3084, 5314, 1620, 1656,
                                                 3504, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5854, 0, 3, 3114, 5359, 1656, 1692,
                                                 3564, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5944, 0, 3, 3144, 5404, 1692, 1728,
                                                 3624, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6034, 0, 3, 3174, 5449, 1728, 1764,
                                                 3684, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6124, 0, 3, 3204, 5494, 1764, 1800,
                                                 3744, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6214, 0, 3, 3264, 5539, 1872, 1908,
                                                 3864, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6304, 0, 3, 3294, 5584, 1908, 1944,
                                                 3924, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6394, 0, 3, 3324, 5629, 1944, 1980,
                                                 3984, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6484, 0, 3, 3354, 5674, 1980, 2016,
                                                 4044, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6574, 0, 3, 3384, 5719, 2016, 2052,
                                                 4104, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6664, 0, 3, 3444, 5764, 2124, 2184,
                                                 4164, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6814, 0, 3, 3504, 5854, 2184, 2244,
                                                 4264, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6964, 0, 3, 3564, 5944, 2244, 2304,
                                                 4364, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7114, 0, 3, 3624, 6034, 2304, 2364,
                                                 4464, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7264, 0, 3, 3684, 6124, 2364, 2424,
                                                 4564, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7414, 0, 3, 3804, 6214, 2544, 2604,
                                                 4664, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7564, 0, 3, 3864, 6304, 2604, 2664,
                                                 4764, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7714, 0, 3, 3924, 6394, 2664, 2724,
                                                 4864, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7864, 0, 3, 3984, 6484, 2724, 2784,
                                                 4964, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8014, 0, 3, 4044, 6574, 2784, 2844,
                                                 5064, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8164, 3, 2964, 2974, 5179, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8185, 3, 2974, 2984, 5194, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8206, 3, 2984, 2994, 5209, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8227, 3, 2994, 3004, 5224, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8248, 3, 3024, 3034, 5254, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8269, 3, 3034, 3044, 5269, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8290, 3, 3044, 3054, 5284, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8311, 3, 3054, 3064, 5299, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 8332, 0, 3, 5164, 8164, 5359, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8395, 0, 3, 5179, 8185, 5404, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8458, 0, 3, 5194, 8206, 5449, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8521, 0, 3, 5209, 8227, 5494, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8584, 0, 3, 5239, 8248, 5584, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8647, 0, 3, 5254, 8269, 5629, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8710, 0, 3, 5269, 8290, 5674, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8773, 0, 3, 5284, 8311, 5719, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 8836, 0, 3, 5314, 8332, 3444, 3504,
                                                 5854, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8962, 0, 3, 5359, 8395, 3504, 3564,
                                                 5944, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9088, 0, 3, 5404, 8458, 3564, 3624,
                                                 6034, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9214, 0, 3, 5449, 8521, 3624, 3684,
                                                 6124, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9340, 0, 3, 5539, 8584, 3804, 3864,
                                                 6304, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9466, 0, 3, 5584, 8647, 3864, 3924,
                                                 6394, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9592, 0, 3, 5629, 8710, 3924, 3984,
                                                 6484, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9718, 0, 3, 5674, 8773, 3984, 4044,
                                                 6574, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9844, 0, 3, 5854, 8962, 4164, 4264,
                                                 6964, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10054, 0, 3, 5944, 9088, 4264, 4364,
                                                 7114, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10264, 0, 3, 6034, 9214, 4364, 4464,
                                                 7264, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10474, 0, 3, 6304, 9466, 4664, 4764,
                                                 7714, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10684, 0, 3, 6394, 9592, 4764, 4864,
                                                 7864, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10894, 0, 3, 6484, 9718, 4864, 4964,
                                                 8014, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11104, 3, 5164, 5179, 8185, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11132, 3, 5179, 5194, 8206, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11160, 3, 5194, 5209, 8227, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11188, 3, 5239, 5254, 8269, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11216, 3, 5254, 5269, 8290, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11244, 3, 5269, 5284, 8311, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 11272, 0, 3, 8164, 11104, 8395, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11356, 0, 3, 8185, 11132, 8458, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11440, 0, 3, 8206, 11160, 8521, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11524, 0, 3, 8248, 11188, 8647, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11608, 0, 3, 8269, 11216, 8710, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11692, 0, 3, 8290, 11244, 8773, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 11776, 0, 3, 8332, 11272, 5764, 5854,
                                                 8962, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11944, 0, 3, 8395, 11356, 5854, 5944,
                                                 9088, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12112, 0, 3, 8458, 11440, 5944, 6034,
                                                 9214, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12280, 0, 3, 8584, 11524, 6214, 6304,
                                                 9466, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12448, 0, 3, 8647, 11608, 6304, 6394,
                                                 9592, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12616, 0, 3, 8710, 11692, 6394, 6484,
                                                 9718, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12784, 0, 3, 8836, 11776, 6664, 6814,
                                                 9844, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13064, 0, 3, 8962, 11944, 6814, 6964,
                                                 10054, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13344, 0, 3, 9088, 12112, 6964, 7114,
                                                 10264, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13624, 0, 3, 9340, 12280, 7414, 7564,
                                                 10474, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 13904, 0, 3, 9466, 12448, 7564, 7714,
                                                 10684, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 14184, 0, 3, 9592, 12616, 7714, 7864,
                                                 10894, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14464, 3, 8164, 8185, 11132, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14500, 3, 8185, 8206, 11160, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14536, 3, 8248, 8269, 11216, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14572, 3, 8269, 8290, 11244, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 14608, 0, 3, 11104, 14464, 11356, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14716, 0, 3, 11132, 14500, 11440, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14824, 0, 3, 11188, 14536, 11608, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14932, 0, 3, 11216, 14572, 11692, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 15040, 0, 3, 11272, 14608, 8836, 8962,
                                                 11944, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 15256, 0, 3, 11356, 14716, 8962, 9088,
                                                 12112, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 15472, 0, 3, 11524, 14824, 9340, 9466,
                                                 12448, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 15688, 0, 3, 11608, 14932, 9466, 9592,
                                                 12616, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 15904, 0, 3, 11944, 15256, 9844, 10054,
                                                 13344, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 16264, 0, 3, 12448, 15688, 10474, 10684,
                                                 14184, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 16624, 3, 11104, 11132, 14500, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 16669, 3, 11188, 11216, 14572, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 16714, 0, 3, 14464, 16624, 14716, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 16849, 0, 3, 14536, 16669, 14932, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 16984, 0, 3, 14608, 16714, 11776, 11944,
                                                 15256, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 17254, 0, 3, 14824, 16849, 12280, 12448,
                                                 15688, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 17524, 0, 3, 15040, 16984, 12784, 13064,
                                                 15904, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 17974, 0, 3, 15472, 17254, 13624, 13904,
                                                 16264, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 18424, 17524, 900, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 19324, 18874, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 19324, 17, nmax);

    simdtrf::transform_l_inner(buffer, 19324, 18424, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 119 * nvalues, nvalues, buffer, 19324, 17, nmax);
}

}  // namespace simdt2ceri
