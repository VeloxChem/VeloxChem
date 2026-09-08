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


#include "SimdElectronRepulsionRecHL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGL.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHL.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hl_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(33467, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 9, 10, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 10, 11, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 11, 12, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 12, 13, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 13, 14, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 14, 15, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 15, 16, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 16, 17, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 119, 0, 17, 18, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 20, 23, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 23, 26, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 26, 29, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 29, 32, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 32, 35, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 35, 38, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 38, 41, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 41, 44, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 44, 47, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 47, 50, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 50, 53, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 235, 0, 59, 65, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 250, 0, 65, 71, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 265, 0, 71, 77, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 280, 0, 77, 83, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 295, 0, 83, 89, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 310, 0, 89, 95, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 325, 0, 95, 101, 205, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 340, 0, 101, 107, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 355, 0, 107, 113, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 370, 0, 125, 135, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 391, 0, 135, 145, 250, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 412, 0, 145, 155, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 433, 0, 155, 165, 280, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 454, 0, 165, 175, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 475, 0, 175, 185, 310, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 496, 0, 185, 195, 325, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 517, 0, 195, 205, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 538, 0, 205, 215, 355, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 559, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 562, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 565, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 568, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 571, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 574, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 577, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 580, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 583, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 586, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 589, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 598, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 607, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 616, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 625, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 634, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 643, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 652, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 661, 3, 17, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 670, 3, 18, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 679, 0, 3, 26, 589, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 697, 0, 3, 29, 598, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 715, 0, 3, 32, 607, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 733, 0, 3, 35, 616, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 751, 0, 3, 38, 625, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 769, 0, 3, 41, 634, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 787, 0, 3, 44, 643, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 805, 0, 3, 47, 652, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 823, 0, 3, 50, 661, 113, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 841, 0, 3, 53, 670, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 859, 0, 3, 65, 697, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 889, 0, 3, 71, 715, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 919, 0, 3, 77, 733, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 949, 0, 3, 83, 751, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 979, 0, 3, 89, 769, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1009, 0, 3, 95, 787, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1039, 0, 3, 101, 805, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1069, 0, 3, 107, 823, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1099, 0, 3, 113, 841, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1129, 0, 3, 145, 889, 250, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1174, 0, 3, 155, 919, 265, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1219, 0, 3, 165, 949, 280, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1264, 0, 3, 175, 979, 295, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1309, 0, 3, 185, 1009, 310, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1354, 0, 3, 195, 1039, 325, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1399, 0, 3, 205, 1069, 340, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1444, 0, 3, 215, 1099, 355, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1489, 0, 3, 250, 1174, 412, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1552, 0, 3, 265, 1219, 433, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1615, 0, 3, 280, 1264, 454, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1678, 0, 3, 295, 1309, 475, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1741, 0, 3, 310, 1354, 496, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1804, 0, 3, 325, 1399, 517, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1867, 0, 3, 340, 1444, 538, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1930, 3, 9, 10, 562, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1936, 3, 10, 11, 565, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1942, 3, 11, 12, 568, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1948, 3, 12, 13, 571, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1954, 3, 13, 14, 574, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1960, 3, 14, 15, 577, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1966, 3, 15, 16, 580, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1972, 3, 16, 17, 583, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1978, 3, 17, 18, 586, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1984, 0, 3, 559, 1930, 598, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2002, 0, 3, 562, 1936, 607, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2020, 0, 3, 565, 1942, 616, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2038, 0, 3, 568, 1948, 625, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2056, 0, 3, 571, 1954, 634, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2074, 0, 3, 574, 1960, 643, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2092, 0, 3, 577, 1966, 652, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2110, 0, 3, 580, 1972, 661, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2128, 0, 3, 583, 1978, 670, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2146, 0, 3, 589, 1984, 59, 65, 697,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2182, 0, 3, 598, 2002, 65, 71, 715,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2218, 0, 3, 607, 2020, 71, 77, 733,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2254, 0, 3, 616, 2038, 77, 83, 751,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2290, 0, 3, 625, 2056, 83, 89, 769,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2326, 0, 3, 634, 2074, 89, 95, 787,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2362, 0, 3, 643, 2092, 95, 101, 805,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2398, 0, 3, 652, 2110, 101, 107, 823,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2434, 0, 3, 661, 2128, 107, 113, 841,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2470, 0, 3, 679, 2146, 125, 135, 859,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2530, 0, 3, 697, 2182, 135, 145, 889,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2590, 0, 3, 715, 2218, 145, 155, 919,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2650, 0, 3, 733, 2254, 155, 165, 949,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2710, 0, 3, 751, 2290, 165, 175, 979,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2770, 0, 3, 769, 2326, 175, 185, 1009,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2830, 0, 3, 787, 2362, 185, 195, 1039,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2890, 0, 3, 805, 2398, 195, 205, 1069,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2950, 0, 3, 823, 2434, 205, 215, 1099,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3010, 0, 3, 2146, 2182, 889, 2590, 235,
                                                 250, 1174, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3100, 0, 3, 2182, 2218, 919, 2650, 250,
                                                 265, 1219, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3190, 0, 3, 2218, 2254, 949, 2710, 265,
                                                 280, 1264, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3280, 0, 3, 2254, 2290, 979, 2770, 280,
                                                 295, 1309, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3370, 0, 3, 2290, 2326, 1009, 2830, 295,
                                                 310, 1354, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3460, 0, 3, 2326, 2362, 1039, 2890, 310,
                                                 325, 1399, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3550, 0, 3, 2362, 2398, 1069, 2950, 325,
                                                 340, 1444, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3640, 0, 3, 2470, 2530, 1129, 3010, 370,
                                                 391, 1489, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3766, 0, 3, 2530, 2590, 1174, 3100, 391,
                                                 412, 1552, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3892, 0, 3, 2590, 2650, 1219, 3190, 412,
                                                 433, 1615, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4018, 0, 3, 2650, 2710, 1264, 3280, 433,
                                                 454, 1678, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4144, 0, 3, 2710, 2770, 1309, 3370, 454,
                                                 475, 1741, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4270, 0, 3, 2770, 2830, 1354, 3460, 475,
                                                 496, 1804, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4396, 0, 3, 2830, 2890, 1399, 3550, 496,
                                                 517, 1867, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4522, 3, 559, 562, 1936, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4532, 3, 562, 565, 1942, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4542, 3, 565, 568, 1948, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4552, 3, 568, 571, 1954, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4562, 3, 571, 574, 1960, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4572, 3, 574, 577, 1966, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4582, 3, 577, 580, 1972, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4592, 3, 580, 583, 1978, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4602, 0, 3, 1930, 4522, 2002, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4632, 0, 3, 1936, 4532, 2020, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4662, 0, 3, 1942, 4542, 2038, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4692, 0, 3, 1948, 4552, 2056, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4722, 0, 3, 1954, 4562, 2074, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4752, 0, 3, 1960, 4572, 2092, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4782, 0, 3, 1966, 4582, 2110, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4812, 0, 3, 1972, 4592, 2128, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4842, 0, 3, 1984, 4602, 679, 697, 2182,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4902, 0, 3, 2002, 4632, 697, 715, 2218,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4962, 0, 3, 2020, 4662, 715, 733, 2254,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5022, 0, 3, 2038, 4692, 733, 751, 2290,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5082, 0, 3, 2056, 4722, 751, 769, 2326,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5142, 0, 3, 2074, 4752, 769, 787, 2362,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5202, 0, 3, 2092, 4782, 787, 805, 2398,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5262, 0, 3, 2110, 4812, 805, 823, 2434,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5322, 0, 3, 2182, 4902, 859, 889, 2590,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5422, 0, 3, 2218, 4962, 889, 919, 2650,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5522, 0, 3, 2254, 5022, 919, 949, 2710,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5622, 0, 3, 2290, 5082, 949, 979, 2770,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5722, 0, 3, 2326, 5142, 979, 1009, 2830,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5822, 0, 3, 2362, 5202, 1009, 1039,
                                                 2890, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5922, 0, 3, 2398, 5262, 1039, 1069,
                                                 2950, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6022, 0, 3, 4842, 4902, 2590, 5422,
                                                 1129, 1174, 3100, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6172, 0, 3, 4902, 4962, 2650, 5522,
                                                 1174, 1219, 3190, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6322, 0, 3, 4962, 5022, 2710, 5622,
                                                 1219, 1264, 3280, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6472, 0, 3, 5022, 5082, 2770, 5722,
                                                 1264, 1309, 3370, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6622, 0, 3, 5082, 5142, 2830, 5822,
                                                 1309, 1354, 3460, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6772, 0, 3, 5142, 5202, 2890, 5922,
                                                 1354, 1399, 3550, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6922, 0, 3, 5322, 5422, 3100, 6172,
                                                 1489, 1552, 3892, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7132, 0, 3, 5422, 5522, 3190, 6322,
                                                 1552, 1615, 4018, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7342, 0, 3, 5522, 5622, 3280, 6472,
                                                 1615, 1678, 4144, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7552, 0, 3, 5622, 5722, 3370, 6622,
                                                 1678, 1741, 4270, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7762, 0, 3, 5722, 5822, 3460, 6772,
                                                 1741, 1804, 4396, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7972, 3, 1930, 1936, 4532, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7987, 3, 1936, 1942, 4542, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8002, 3, 1942, 1948, 4552, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8017, 3, 1948, 1954, 4562, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8032, 3, 1954, 1960, 4572, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8047, 3, 1960, 1966, 4582, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 8062, 3, 1966, 1972, 4592, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 8077, 0, 3, 4522, 7972, 4632, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8122, 0, 3, 4532, 7987, 4662, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8167, 0, 3, 4542, 8002, 4692, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8212, 0, 3, 4552, 8017, 4722, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8257, 0, 3, 4562, 8032, 4752, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8302, 0, 3, 4572, 8047, 4782, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 8347, 0, 3, 4582, 8062, 4812, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 8392, 0, 3, 4602, 8077, 2146, 2182,
                                                 4902, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8482, 0, 3, 4632, 8122, 2182, 2218,
                                                 4962, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8572, 0, 3, 4662, 8167, 2218, 2254,
                                                 5022, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8662, 0, 3, 4692, 8212, 2254, 2290,
                                                 5082, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8752, 0, 3, 4722, 8257, 2290, 2326,
                                                 5142, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8842, 0, 3, 4752, 8302, 2326, 2362,
                                                 5202, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8932, 0, 3, 4782, 8347, 2362, 2398,
                                                 5262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9022, 0, 3, 4842, 8392, 2470, 2530,
                                                 5322, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9172, 0, 3, 4902, 8482, 2530, 2590,
                                                 5422, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9322, 0, 3, 4962, 8572, 2590, 2650,
                                                 5522, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9472, 0, 3, 5022, 8662, 2650, 2710,
                                                 5622, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9622, 0, 3, 5082, 8752, 2710, 2770,
                                                 5722, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9772, 0, 3, 5142, 8842, 2770, 2830,
                                                 5822, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9922, 0, 3, 5202, 8932, 2830, 2890,
                                                 5922, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10072, 0, 3, 8392, 8482, 5422, 9322,
                                                 3010, 3100, 6172, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10297, 0, 3, 8482, 8572, 5522, 9472,
                                                 3100, 3190, 6322, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10522, 0, 3, 8572, 8662, 5622, 9622,
                                                 3190, 3280, 6472, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10747, 0, 3, 8662, 8752, 5722, 9772,
                                                 3280, 3370, 6622, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10972, 0, 3, 8752, 8842, 5822, 9922,
                                                 3370, 3460, 6772, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 11197, 0, 3, 9022, 9172, 6022, 10072,
                                                 3640, 3766, 6922, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 11512, 0, 3, 9172, 9322, 6172, 10297,
                                                 3766, 3892, 7132, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 11827, 0, 3, 9322, 9472, 6322, 10522,
                                                 3892, 4018, 7342, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12142, 0, 3, 9472, 9622, 6472, 10747,
                                                 4018, 4144, 7552, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12457, 0, 3, 9622, 9772, 6622, 10972,
                                                 4144, 4270, 7762, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12772, 3, 4522, 4532, 7987, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12793, 3, 4532, 4542, 8002, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12814, 3, 4542, 4552, 8017, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12835, 3, 4552, 4562, 8032, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12856, 3, 4562, 4572, 8047, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 12877, 3, 4572, 4582, 8062, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 12898, 0, 3, 7972, 12772, 8122, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 12961, 0, 3, 7987, 12793, 8167, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13024, 0, 3, 8002, 12814, 8212, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13087, 0, 3, 8017, 12835, 8257, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13150, 0, 3, 8032, 12856, 8302, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13213, 0, 3, 8047, 12877, 8347, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 13276, 0, 3, 8077, 12898, 4842, 4902,
                                                 8482, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13402, 0, 3, 8122, 12961, 4902, 4962,
                                                 8572, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13528, 0, 3, 8167, 13024, 4962, 5022,
                                                 8662, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13654, 0, 3, 8212, 13087, 5022, 5082,
                                                 8752, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13780, 0, 3, 8257, 13150, 5082, 5142,
                                                 8842, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13906, 0, 3, 8302, 13213, 5142, 5202,
                                                 8932, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14032, 0, 3, 8482, 13402, 5322, 5422,
                                                 9322, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14242, 0, 3, 8572, 13528, 5422, 5522,
                                                 9472, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14452, 0, 3, 8662, 13654, 5522, 5622,
                                                 9622, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14662, 0, 3, 8752, 13780, 5622, 5722,
                                                 9772, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14872, 0, 3, 8842, 13906, 5722, 5822,
                                                 9922, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 15082, 0, 3, 13276, 13402, 9322, 14242,
                                                 6022, 6172, 10297, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 15397, 0, 3, 13402, 13528, 9472, 14452,
                                                 6172, 6322, 10522, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 15712, 0, 3, 13528, 13654, 9622, 14662,
                                                 6322, 6472, 10747, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16027, 0, 3, 13654, 13780, 9772, 14872,
                                                 6472, 6622, 10972, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 16342, 0, 3, 14032, 14242, 10297, 15397,
                                                 6922, 7132, 11827, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 16783, 0, 3, 14242, 14452, 10522, 15712,
                                                 7132, 7342, 12142, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 17224, 0, 3, 14452, 14662, 10747, 16027,
                                                 7342, 7552, 12457, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17665, 3, 7972, 7987, 12793, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17693, 3, 7987, 8002, 12814, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17721, 3, 8002, 8017, 12835, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17749, 3, 8017, 8032, 12856, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 17777, 3, 8032, 8047, 12877, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 17805, 0, 3, 12772, 17665, 12961, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 17889, 0, 3, 12793, 17693, 13024, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 17973, 0, 3, 12814, 17721, 13087, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18057, 0, 3, 12835, 17749, 13150, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 18141, 0, 3, 12856, 17777, 13213, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 18225, 0, 3, 12898, 17805, 8392, 8482,
                                                 13402, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 18393, 0, 3, 12961, 17889, 8482, 8572,
                                                 13528, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 18561, 0, 3, 13024, 17973, 8572, 8662,
                                                 13654, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 18729, 0, 3, 13087, 18057, 8662, 8752,
                                                 13780, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 18897, 0, 3, 13150, 18141, 8752, 8842,
                                                 13906, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 19065, 0, 3, 13276, 18225, 9022, 9172,
                                                 14032, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 19345, 0, 3, 13402, 18393, 9172, 9322,
                                                 14242, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 19625, 0, 3, 13528, 18561, 9322, 9472,
                                                 14452, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 19905, 0, 3, 13654, 18729, 9472, 9622,
                                                 14662, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 20185, 0, 3, 13780, 18897, 9622, 9772,
                                                 14872, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 20465, 0, 3, 18225, 18393, 14242, 19625,
                                                 10072, 10297, 15397, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 20885, 0, 3, 18393, 18561, 14452, 19905,
                                                 10297, 10522, 15712, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 21305, 0, 3, 18561, 18729, 14662, 20185,
                                                 10522, 10747, 16027, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 21725, 0, 3, 19065, 19345, 15082, 20465,
                                                 11197, 11512, 16342, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 22313, 0, 3, 19345, 19625, 15397, 20885,
                                                 11512, 11827, 16783, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 22901, 0, 3, 19625, 19905, 15712, 21305,
                                                 11827, 12142, 17224, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 23489, 3, 12772, 12793, 17693, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 23525, 3, 12793, 12814, 17721, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 23561, 3, 12814, 12835, 17749, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 23597, 3, 12835, 12856, 17777, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 23633, 0, 3, 17665, 23489, 17889, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 23741, 0, 3, 17693, 23525, 17973, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 23849, 0, 3, 17721, 23561, 18057, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 23957, 0, 3, 17749, 23597, 18141, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 24065, 0, 3, 17805, 23633, 13276, 13402,
                                                 18393, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 24281, 0, 3, 17889, 23741, 13402, 13528,
                                                 18561, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 24497, 0, 3, 17973, 23849, 13528, 13654,
                                                 18729, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 24713, 0, 3, 18057, 23957, 13654, 13780,
                                                 18897, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 24929, 0, 3, 18393, 24281, 14032, 14242,
                                                 19625, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 25289, 0, 3, 18561, 24497, 14242, 14452,
                                                 19905, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 25649, 0, 3, 18729, 24713, 14452, 14662,
                                                 20185, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 26009, 0, 3, 24065, 24281, 19625, 25289,
                                                 15082, 15397, 20885, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 26549, 0, 3, 24281, 24497, 19905, 25649,
                                                 15397, 15712, 21305, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 27089, 0, 3, 24929, 25289, 20885, 26549,
                                                 16342, 16783, 22901, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 27845, 3, 17665, 17693, 23525, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 27890, 3, 17693, 17721, 23561, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 27935, 3, 17721, 17749, 23597, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 27980, 0, 3, 23489, 27845, 23741, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 28115, 0, 3, 23525, 27890, 23849, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 28250, 0, 3, 23561, 27935, 23957, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 28385, 0, 3, 23633, 27980, 18225, 18393,
                                                 24281, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 28655, 0, 3, 23741, 28115, 18393, 18561,
                                                 24497, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 28925, 0, 3, 23849, 28250, 18561, 18729,
                                                 24713, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 29195, 0, 3, 24065, 28385, 19065, 19345,
                                                 24929, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 29645, 0, 3, 24281, 28655, 19345, 19625,
                                                 25289, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 30095, 0, 3, 24497, 28925, 19625, 19905,
                                                 25649, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 30545, 0, 3, 28385, 28655, 25289, 30095,
                                                 20465, 20885, 26549, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 31220, 0, 3, 29195, 29645, 26009, 30545,
                                                 21725, 22313, 27089, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 32165, 31220, 945, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 33110, 32165, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 33110, 17, nmax);
}

}  // namespace simdt2ceri
