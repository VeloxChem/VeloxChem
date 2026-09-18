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


#include "SimdElectronRepulsionRsRecFK.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_fk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 13088, 12218, 720, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 318, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 321, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 324, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 327, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 330, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 333, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 336, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 339, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 342, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 345, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 348, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 351, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 354, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 357, 3, 27, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 360, 3, 9, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 369, 3, 10, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 378, 3, 11, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 387, 3, 12, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 396, 3, 13, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 405, 3, 14, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 414, 3, 15, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 423, 3, 20, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 432, 3, 21, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 441, 3, 22, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 450, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 459, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 468, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 477, 3, 26, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 486, 0, 3, 31, 360, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 504, 0, 3, 34, 369, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 522, 0, 3, 37, 378, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 540, 0, 3, 40, 387, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 558, 0, 3, 43, 396, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 576, 0, 3, 46, 405, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 594, 0, 3, 49, 414, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 612, 0, 3, 58, 423, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 630, 0, 3, 61, 432, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 648, 0, 3, 64, 441, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 666, 0, 3, 67, 450, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 684, 0, 3, 70, 459, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 702, 0, 3, 73, 468, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 720, 0, 3, 76, 477, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 738, 0, 3, 82, 486, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 768, 0, 3, 88, 504, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 798, 0, 3, 94, 522, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 828, 0, 3, 100, 540, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 858, 0, 3, 106, 558, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 888, 0, 3, 112, 576, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 918, 0, 3, 118, 594, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 948, 0, 3, 130, 612, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 978, 0, 3, 136, 630, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1008, 0, 3, 142, 648, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1038, 0, 3, 148, 666, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1068, 0, 3, 154, 684, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1098, 0, 3, 160, 702, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1128, 0, 3, 166, 720, 308, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1158, 3, 9, 10, 321, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1164, 3, 10, 11, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1170, 3, 11, 12, 327, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1176, 3, 12, 13, 330, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1182, 3, 13, 14, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1188, 3, 14, 15, 336, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1194, 3, 20, 21, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1200, 3, 21, 22, 345, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1206, 3, 22, 23, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1212, 3, 23, 24, 351, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1218, 3, 24, 25, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1224, 3, 25, 26, 357, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1230, 0, 3, 318, 1158, 369, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1248, 0, 3, 321, 1164, 378, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1266, 0, 3, 324, 1170, 387, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1284, 0, 3, 327, 1176, 396, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1302, 0, 3, 330, 1182, 405, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1320, 0, 3, 333, 1188, 414, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1338, 0, 3, 339, 1194, 432, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1356, 0, 3, 342, 1200, 441, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1374, 0, 3, 345, 1206, 450, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1392, 0, 3, 348, 1212, 459, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1410, 0, 3, 351, 1218, 468, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1428, 0, 3, 354, 1224, 477, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1446, 0, 3, 360, 1230, 82, 88, 504,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1482, 0, 3, 369, 1248, 88, 94, 522,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1518, 0, 3, 378, 1266, 94, 100, 540,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1554, 0, 3, 387, 1284, 100, 106, 558,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1590, 0, 3, 396, 1302, 106, 112, 576,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1626, 0, 3, 405, 1320, 112, 118, 594,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1662, 0, 3, 423, 1338, 130, 136, 630,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1698, 0, 3, 432, 1356, 136, 142, 648,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1734, 0, 3, 441, 1374, 142, 148, 666,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1770, 0, 3, 450, 1392, 148, 154, 684,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1806, 0, 3, 459, 1410, 154, 160, 702,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1842, 0, 3, 468, 1428, 160, 166, 720,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1878, 0, 3, 504, 1482, 178, 188, 798,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1938, 0, 3, 522, 1518, 188, 198, 828,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1998, 0, 3, 540, 1554, 198, 208, 858,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2058, 0, 3, 558, 1590, 208, 218, 888,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2118, 0, 3, 576, 1626, 218, 228, 918,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2178, 0, 3, 630, 1698, 248, 258, 1008,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2238, 0, 3, 648, 1734, 258, 268, 1038,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2298, 0, 3, 666, 1770, 268, 278, 1068,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2358, 0, 3, 684, 1806, 278, 288, 1098,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2418, 0, 3, 702, 1842, 288, 298, 1128,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2478, 3, 318, 321, 1164, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2488, 3, 321, 324, 1170, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2498, 3, 324, 327, 1176, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2508, 3, 327, 330, 1182, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2518, 3, 330, 333, 1188, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2528, 3, 339, 342, 1200, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2538, 3, 342, 345, 1206, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2548, 3, 345, 348, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2558, 3, 348, 351, 1218, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2568, 3, 351, 354, 1224, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2578, 0, 3, 1158, 2478, 1248, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2608, 0, 3, 1164, 2488, 1266, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2638, 0, 3, 1170, 2498, 1284, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2668, 0, 3, 1176, 2508, 1302, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2698, 0, 3, 1182, 2518, 1320, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2728, 0, 3, 1194, 2528, 1356, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2758, 0, 3, 1200, 2538, 1374, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2788, 0, 3, 1206, 2548, 1392, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2818, 0, 3, 1212, 2558, 1410, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2848, 0, 3, 1218, 2568, 1428, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2878, 0, 3, 1230, 2578, 486, 504, 1482,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2938, 0, 3, 1248, 2608, 504, 522, 1518,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2998, 0, 3, 1266, 2638, 522, 540, 1554,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3058, 0, 3, 1284, 2668, 540, 558, 1590,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3118, 0, 3, 1302, 2698, 558, 576, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3178, 0, 3, 1338, 2728, 612, 630, 1698,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3238, 0, 3, 1356, 2758, 630, 648, 1734,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3298, 0, 3, 1374, 2788, 648, 666, 1770,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3358, 0, 3, 1392, 2818, 666, 684, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3418, 0, 3, 1410, 2848, 684, 702, 1842,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3478, 0, 3, 1446, 2878, 738, 768, 1878,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3578, 0, 3, 1482, 2938, 768, 798, 1938,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3678, 0, 3, 1518, 2998, 798, 828, 1998,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3778, 0, 3, 1554, 3058, 828, 858, 2058,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3878, 0, 3, 1590, 3118, 858, 888, 2118,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3978, 0, 3, 1662, 3178, 948, 978, 2178,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4078, 0, 3, 1698, 3238, 978, 1008, 2238,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4178, 0, 3, 1734, 3298, 1008, 1038,
                                                 2298, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4278, 0, 3, 1770, 3358, 1038, 1068,
                                                 2358, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4378, 0, 3, 1806, 3418, 1068, 1098,
                                                 2418, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4478, 3, 1158, 1164, 2488, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4493, 3, 1164, 1170, 2498, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4508, 3, 1170, 1176, 2508, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4523, 3, 1176, 1182, 2518, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4538, 3, 1194, 1200, 2538, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4553, 3, 1200, 1206, 2548, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4568, 3, 1206, 1212, 2558, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4583, 3, 1212, 1218, 2568, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4598, 0, 3, 2478, 4478, 2608, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4643, 0, 3, 2488, 4493, 2638, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4688, 0, 3, 2498, 4508, 2668, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4733, 0, 3, 2508, 4523, 2698, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4778, 0, 3, 2528, 4538, 2758, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4823, 0, 3, 2538, 4553, 2788, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4868, 0, 3, 2548, 4568, 2818, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4913, 0, 3, 2558, 4583, 2848, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 4958, 0, 3, 2578, 4598, 1446, 1482,
                                                 2938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5048, 0, 3, 2608, 4643, 1482, 1518,
                                                 2998, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5138, 0, 3, 2638, 4688, 1518, 1554,
                                                 3058, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5228, 0, 3, 2668, 4733, 1554, 1590,
                                                 3118, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5318, 0, 3, 2728, 4778, 1662, 1698,
                                                 3238, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5408, 0, 3, 2758, 4823, 1698, 1734,
                                                 3298, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5498, 0, 3, 2788, 4868, 1734, 1770,
                                                 3358, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5588, 0, 3, 2818, 4913, 1770, 1806,
                                                 3418, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5678, 0, 3, 2938, 5048, 1878, 1938,
                                                 3678, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5828, 0, 3, 2998, 5138, 1938, 1998,
                                                 3778, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5978, 0, 3, 3058, 5228, 1998, 2058,
                                                 3878, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6128, 0, 3, 3238, 5408, 2178, 2238,
                                                 4178, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6278, 0, 3, 3298, 5498, 2238, 2298,
                                                 4278, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6428, 0, 3, 3358, 5588, 2298, 2358,
                                                 4378, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6578, 3, 2478, 2488, 4493, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6599, 3, 2488, 2498, 4508, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6620, 3, 2498, 2508, 4523, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6641, 3, 2528, 2538, 4553, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6662, 3, 2538, 2548, 4568, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6683, 3, 2548, 2558, 4583, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 6704, 0, 3, 4478, 6578, 4643, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6767, 0, 3, 4493, 6599, 4688, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6830, 0, 3, 4508, 6620, 4733, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6893, 0, 3, 4538, 6641, 4823, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6956, 0, 3, 4553, 6662, 4868, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7019, 0, 3, 4568, 6683, 4913, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 7082, 0, 3, 4598, 6704, 2878, 2938,
                                                 5048, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7208, 0, 3, 4643, 6767, 2938, 2998,
                                                 5138, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7334, 0, 3, 4688, 6830, 2998, 3058,
                                                 5228, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7460, 0, 3, 4778, 6893, 3178, 3238,
                                                 5408, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7586, 0, 3, 4823, 6956, 3238, 3298,
                                                 5498, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7712, 0, 3, 4868, 7019, 3298, 3358,
                                                 5588, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7838, 0, 3, 4958, 7082, 3478, 3578,
                                                 5678, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8048, 0, 3, 5048, 7208, 3578, 3678,
                                                 5828, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8258, 0, 3, 5138, 7334, 3678, 3778,
                                                 5978, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8468, 0, 3, 5318, 7460, 3978, 4078,
                                                 6128, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8678, 0, 3, 5408, 7586, 4078, 4178,
                                                 6278, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8888, 0, 3, 5498, 7712, 4178, 4278,
                                                 6428, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9098, 3, 4478, 4493, 6599, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9126, 3, 4493, 4508, 6620, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9154, 3, 4538, 4553, 6662, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9182, 3, 4553, 4568, 6683, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 9210, 0, 3, 6578, 9098, 6767, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9294, 0, 3, 6599, 9126, 6830, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9378, 0, 3, 6641, 9154, 6956, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9462, 0, 3, 6662, 9182, 7019, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 9546, 0, 3, 6704, 9210, 4958, 5048,
                                                 7208, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 9714, 0, 3, 6767, 9294, 5048, 5138,
                                                 7334, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 9882, 0, 3, 6893, 9378, 5318, 5408,
                                                 7586, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 10050, 0, 3, 6956, 9462, 5408, 5498,
                                                 7712, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 10218, 0, 3, 7208, 9714, 5678, 5828,
                                                 8258, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 10498, 0, 3, 7586, 10050, 6128, 6278,
                                                 8888, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 10778, 3, 6578, 6599, 9126, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 10814, 3, 6641, 6662, 9182, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 10850, 0, 3, 9098, 10778, 9294, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 10958, 0, 3, 9154, 10814, 9462, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 11066, 0, 3, 9210, 10850, 7082, 7208,
                                                 9714, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 11282, 0, 3, 9378, 10958, 7460, 7586,
                                                 10050, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 11498, 0, 3, 9546, 11066, 7838, 8048,
                                                 10218, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 11858, 0, 3, 9882, 11282, 8468, 8678,
                                                 10498, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 12218, 11498, 720, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 12938, 12578, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 12938, 15, nmax);

    simdtrf::transform_k_inner(buffer, 12938, 12218, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 105 * nvalues, nvalues, buffer, 12938, 15, nmax);
}

}  // namespace simdt2ceri
