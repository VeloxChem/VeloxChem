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


#include "SimdElectronRepulsionRsRecGH.hpp"

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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_gh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 10941, 10146, 630, nvalues);

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

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 50, 53, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 53, 56, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 56, 59, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 59, 62, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 62, 65, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 65, 68, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 278, 0, 74, 80, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 293, 0, 80, 86, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 308, 0, 86, 92, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 323, 0, 92, 98, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 338, 0, 98, 104, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 353, 0, 116, 122, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 368, 0, 122, 128, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 383, 0, 128, 134, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 398, 0, 134, 140, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 413, 0, 140, 146, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 428, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 431, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 434, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 437, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 440, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 443, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 446, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 449, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 452, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 455, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 458, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 461, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 464, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 467, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 470, 3, 8, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 479, 3, 9, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 488, 3, 10, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 497, 3, 11, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 506, 3, 12, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 515, 3, 13, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 524, 3, 14, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 533, 3, 18, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 542, 3, 19, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 551, 3, 20, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 560, 3, 21, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 569, 3, 22, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 578, 3, 23, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 587, 3, 24, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 596, 0, 3, 26, 470, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 614, 0, 3, 29, 479, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 632, 0, 3, 32, 488, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 650, 0, 3, 35, 497, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 668, 0, 3, 38, 506, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 686, 0, 3, 41, 515, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 704, 0, 3, 44, 524, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 722, 0, 3, 50, 533, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 740, 0, 3, 53, 542, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 758, 0, 3, 56, 551, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 776, 0, 3, 59, 560, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 794, 0, 3, 62, 569, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 812, 0, 3, 65, 578, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 830, 0, 3, 68, 587, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 848, 0, 3, 80, 632, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 878, 0, 3, 86, 650, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 908, 0, 3, 92, 668, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 938, 0, 3, 98, 686, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 968, 0, 3, 104, 704, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 998, 0, 3, 122, 758, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1028, 0, 3, 128, 776, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1058, 0, 3, 134, 794, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1088, 0, 3, 140, 812, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1118, 0, 3, 146, 830, 268, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1148, 0, 3, 158, 848, 278, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1193, 0, 3, 168, 878, 293, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1238, 0, 3, 178, 908, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1283, 0, 3, 188, 938, 323, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1328, 0, 3, 198, 968, 338, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1373, 0, 3, 218, 998, 353, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1418, 0, 3, 228, 1028, 368, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1463, 0, 3, 238, 1058, 383, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1508, 0, 3, 248, 1088, 398, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1553, 0, 3, 258, 1118, 413, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1598, 3, 8, 9, 431, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1604, 3, 9, 10, 434, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1610, 3, 10, 11, 437, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1616, 3, 11, 12, 440, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1622, 3, 12, 13, 443, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1628, 3, 13, 14, 446, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1634, 3, 18, 19, 452, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1640, 3, 19, 20, 455, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1646, 3, 20, 21, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1652, 3, 21, 22, 461, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1658, 3, 22, 23, 464, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1664, 3, 23, 24, 467, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1670, 0, 3, 428, 1598, 479, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1688, 0, 3, 431, 1604, 488, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1706, 0, 3, 434, 1610, 497, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1724, 0, 3, 437, 1616, 506, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1742, 0, 3, 440, 1622, 515, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1760, 0, 3, 443, 1628, 524, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1778, 0, 3, 449, 1634, 542, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1796, 0, 3, 452, 1640, 551, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1814, 0, 3, 455, 1646, 560, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1832, 0, 3, 458, 1652, 569, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1850, 0, 3, 461, 1658, 578, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1868, 0, 3, 464, 1664, 587, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1886, 0, 3, 479, 1688, 74, 80, 632,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1922, 0, 3, 488, 1706, 80, 86, 650,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1958, 0, 3, 497, 1724, 86, 92, 668,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1994, 0, 3, 506, 1742, 92, 98, 686,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2030, 0, 3, 515, 1760, 98, 104, 704,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2066, 0, 3, 542, 1796, 116, 122, 758,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2102, 0, 3, 551, 1814, 122, 128, 776,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2138, 0, 3, 560, 1832, 128, 134, 794,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2174, 0, 3, 569, 1850, 134, 140, 812,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2210, 0, 3, 578, 1868, 140, 146, 830,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2246, 0, 3, 632, 1922, 158, 168, 878,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2306, 0, 3, 650, 1958, 168, 178, 908,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2366, 0, 3, 668, 1994, 178, 188, 938,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2426, 0, 3, 686, 2030, 188, 198, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2486, 0, 3, 758, 2102, 218, 228, 1028,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2546, 0, 3, 776, 2138, 228, 238, 1058,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2606, 0, 3, 794, 2174, 238, 248, 1088,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2666, 0, 3, 812, 2210, 248, 258, 1118,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2726, 0, 3, 1886, 1922, 878, 2306, 278,
                                                 293, 1238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2816, 0, 3, 1922, 1958, 908, 2366, 293,
                                                 308, 1283, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2906, 0, 3, 1958, 1994, 938, 2426, 308,
                                                 323, 1328, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2996, 0, 3, 2066, 2102, 1028, 2546, 353,
                                                 368, 1463, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3086, 0, 3, 2102, 2138, 1058, 2606, 368,
                                                 383, 1508, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3176, 0, 3, 2138, 2174, 1088, 2666, 383,
                                                 398, 1553, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3266, 3, 428, 431, 1604, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3276, 3, 431, 434, 1610, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3286, 3, 434, 437, 1616, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3296, 3, 437, 440, 1622, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3306, 3, 440, 443, 1628, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3316, 3, 449, 452, 1640, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3326, 3, 452, 455, 1646, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3336, 3, 455, 458, 1652, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3346, 3, 458, 461, 1658, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3356, 3, 461, 464, 1664, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3366, 0, 3, 1598, 3266, 1688, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3396, 0, 3, 1604, 3276, 1706, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3426, 0, 3, 1610, 3286, 1724, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3456, 0, 3, 1616, 3296, 1742, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3486, 0, 3, 1622, 3306, 1760, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3516, 0, 3, 1634, 3316, 1796, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3546, 0, 3, 1640, 3326, 1814, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3576, 0, 3, 1646, 3336, 1832, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3606, 0, 3, 1652, 3346, 1850, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3636, 0, 3, 1658, 3356, 1868, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3666, 0, 3, 1670, 3366, 596, 614, 1886,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3726, 0, 3, 1688, 3396, 614, 632, 1922,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3786, 0, 3, 1706, 3426, 632, 650, 1958,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3846, 0, 3, 1724, 3456, 650, 668, 1994,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3906, 0, 3, 1742, 3486, 668, 686, 2030,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3966, 0, 3, 1778, 3516, 722, 740, 2066,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4026, 0, 3, 1796, 3546, 740, 758, 2102,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4086, 0, 3, 1814, 3576, 758, 776, 2138,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4146, 0, 3, 1832, 3606, 776, 794, 2174,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4206, 0, 3, 1850, 3636, 794, 812, 2210,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4266, 0, 3, 1922, 3786, 848, 878, 2306,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4366, 0, 3, 1958, 3846, 878, 908, 2366,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4466, 0, 3, 1994, 3906, 908, 938, 2426,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4566, 0, 3, 2102, 4086, 998, 1028, 2546,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4666, 0, 3, 2138, 4146, 1028, 1058,
                                                 2606, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4766, 0, 3, 2174, 4206, 1058, 1088,
                                                 2666, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4866, 0, 3, 3666, 3726, 2246, 4266,
                                                 1148, 1193, 2726, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5016, 0, 3, 3726, 3786, 2306, 4366,
                                                 1193, 1238, 2816, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5166, 0, 3, 3786, 3846, 2366, 4466,
                                                 1238, 1283, 2906, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5316, 0, 3, 3966, 4026, 2486, 4566,
                                                 1373, 1418, 2996, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5466, 0, 3, 4026, 4086, 2546, 4666,
                                                 1418, 1463, 3086, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5616, 0, 3, 4086, 4146, 2606, 4766,
                                                 1463, 1508, 3176, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5766, 3, 1598, 1604, 3276, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5781, 3, 1604, 1610, 3286, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5796, 3, 1610, 1616, 3296, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5811, 3, 1616, 1622, 3306, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5826, 3, 1634, 1640, 3326, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5841, 3, 1640, 1646, 3336, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5856, 3, 1646, 1652, 3346, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5871, 3, 1652, 1658, 3356, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 5886, 0, 3, 3266, 5766, 3396, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5931, 0, 3, 3276, 5781, 3426, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5976, 0, 3, 3286, 5796, 3456, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6021, 0, 3, 3296, 5811, 3486, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6066, 0, 3, 3316, 5826, 3546, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6111, 0, 3, 3326, 5841, 3576, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6156, 0, 3, 3336, 5856, 3606, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 6201, 0, 3, 3346, 5871, 3636, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 6246, 0, 3, 3396, 5931, 1886, 1922,
                                                 3786, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6336, 0, 3, 3426, 5976, 1922, 1958,
                                                 3846, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6426, 0, 3, 3456, 6021, 1958, 1994,
                                                 3906, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6516, 0, 3, 3546, 6111, 2066, 2102,
                                                 4086, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6606, 0, 3, 3576, 6156, 2102, 2138,
                                                 4146, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6696, 0, 3, 3606, 6201, 2138, 2174,
                                                 4206, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6786, 0, 3, 3786, 6336, 2246, 2306,
                                                 4366, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6936, 0, 3, 3846, 6426, 2306, 2366,
                                                 4466, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7086, 0, 3, 4086, 6606, 2486, 2546,
                                                 4666, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7236, 0, 3, 4146, 6696, 2546, 2606,
                                                 4766, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7386, 0, 3, 6246, 6336, 4366, 6936,
                                                 2726, 2816, 5166, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7611, 0, 3, 6516, 6606, 4666, 7236,
                                                 2996, 3086, 5616, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7836, 3, 3266, 3276, 5781, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7857, 3, 3276, 3286, 5796, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7878, 3, 3286, 3296, 5811, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7899, 3, 3316, 3326, 5841, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7920, 3, 3326, 3336, 5856, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7941, 3, 3336, 3346, 5871, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 7962, 0, 3, 5766, 7836, 5931, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8025, 0, 3, 5781, 7857, 5976, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8088, 0, 3, 5796, 7878, 6021, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8151, 0, 3, 5826, 7899, 6111, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8214, 0, 3, 5841, 7920, 6156, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8277, 0, 3, 5856, 7941, 6201, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 8340, 0, 3, 5886, 7962, 3666, 3726,
                                                 6246, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8466, 0, 3, 5931, 8025, 3726, 3786,
                                                 6336, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8592, 0, 3, 5976, 8088, 3786, 3846,
                                                 6426, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8718, 0, 3, 6066, 8151, 3966, 4026,
                                                 6516, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8844, 0, 3, 6111, 8214, 4026, 4086,
                                                 6606, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8970, 0, 3, 6156, 8277, 4086, 4146,
                                                 6696, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9096, 0, 3, 6336, 8592, 4266, 4366,
                                                 6936, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9306, 0, 3, 6606, 8970, 4566, 4666,
                                                 7236, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 9516, 0, 3, 8340, 8466, 6786, 9096,
                                                 4866, 5016, 7386, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 9831, 0, 3, 8718, 8844, 7086, 9306,
                                                 5316, 5466, 7611, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 10146, 9516, 630, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 10776, 10461, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 10776, 11, nmax);

    simdtrf::transform_h_inner(buffer, 10776, 10146, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 99 * nvalues, nvalues, buffer, 10776, 11, nmax);
}

}  // namespace simdt2ceri
