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


#include "SimdElectronRepulsionRecHK.hpp"

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
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 22813, 21742, 756, nvalues);

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
                                            10, 11, 12}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 18, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 14, 15, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 15, 16, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 16, 17, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 19, 22, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 34, 37, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 37, 40, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 40, 43, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 43, 46, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 202, 0, 52, 58, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 217, 0, 58, 64, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 232, 0, 64, 70, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 247, 0, 70, 76, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 262, 0, 76, 82, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 277, 0, 82, 88, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 292, 0, 88, 94, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 307, 0, 94, 100, 192, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 322, 0, 112, 122, 217, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 343, 0, 122, 132, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 364, 0, 132, 142, 247, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 385, 0, 142, 152, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 406, 0, 152, 162, 277, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 427, 0, 162, 172, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 448, 0, 172, 182, 307, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 469, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 472, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 475, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 478, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 481, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 484, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 487, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 490, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 493, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 496, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 505, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 514, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 523, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 532, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 541, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 550, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 559, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 568, 3, 17, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 577, 0, 3, 22, 496, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 595, 0, 3, 25, 505, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 613, 0, 3, 28, 514, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 631, 0, 3, 31, 523, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 649, 0, 3, 34, 532, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 667, 0, 3, 37, 541, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 685, 0, 3, 40, 550, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 703, 0, 3, 43, 559, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 721, 0, 3, 46, 568, 106, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 739, 0, 3, 52, 577, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 769, 0, 3, 58, 595, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 799, 0, 3, 64, 613, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 829, 0, 3, 70, 631, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 859, 0, 3, 76, 649, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 889, 0, 3, 82, 667, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 919, 0, 3, 88, 685, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 949, 0, 3, 94, 703, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 979, 0, 3, 100, 721, 192, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1009, 0, 3, 122, 799, 217, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1054, 0, 3, 132, 829, 232, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1099, 0, 3, 142, 859, 247, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1144, 0, 3, 152, 889, 262, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1189, 0, 3, 162, 919, 277, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1234, 0, 3, 172, 949, 292, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1279, 0, 3, 182, 979, 307, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1324, 0, 3, 202, 1009, 322, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1387, 0, 3, 217, 1054, 343, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1450, 0, 3, 232, 1099, 364, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1513, 0, 3, 247, 1144, 385, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1576, 0, 3, 262, 1189, 406, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1639, 0, 3, 277, 1234, 427, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1702, 0, 3, 292, 1279, 448, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1765, 3, 9, 10, 472, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1771, 3, 10, 11, 475, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1777, 3, 11, 12, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1783, 3, 12, 13, 481, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1789, 3, 13, 14, 484, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1795, 3, 14, 15, 487, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1801, 3, 15, 16, 490, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1807, 3, 16, 17, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1813, 0, 3, 469, 1765, 505, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1831, 0, 3, 472, 1771, 514, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1849, 0, 3, 475, 1777, 523, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1867, 0, 3, 478, 1783, 532, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1885, 0, 3, 481, 1789, 541, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1903, 0, 3, 484, 1795, 550, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1921, 0, 3, 487, 1801, 559, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1939, 0, 3, 490, 1807, 568, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1957, 0, 3, 496, 1813, 52, 58, 595,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1993, 0, 3, 505, 1831, 58, 64, 613,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2029, 0, 3, 514, 1849, 64, 70, 631,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2065, 0, 3, 523, 1867, 70, 76, 649,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2101, 0, 3, 532, 1885, 76, 82, 667,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2137, 0, 3, 541, 1903, 82, 88, 685,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2173, 0, 3, 550, 1921, 88, 94, 703,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2209, 0, 3, 559, 1939, 94, 100, 721,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2245, 0, 3, 595, 1993, 112, 122, 799,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2305, 0, 3, 613, 2029, 122, 132, 829,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2365, 0, 3, 631, 2065, 132, 142, 859,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2425, 0, 3, 649, 2101, 142, 152, 889,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2485, 0, 3, 667, 2137, 152, 162, 919,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2545, 0, 3, 685, 2173, 162, 172, 949,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2605, 0, 3, 703, 2209, 172, 182, 979,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2665, 0, 3, 1957, 1993, 799, 2305, 202,
                                                 217, 1054, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2755, 0, 3, 1993, 2029, 829, 2365, 217,
                                                 232, 1099, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2845, 0, 3, 2029, 2065, 859, 2425, 232,
                                                 247, 1144, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2935, 0, 3, 2065, 2101, 889, 2485, 247,
                                                 262, 1189, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3025, 0, 3, 2101, 2137, 919, 2545, 262,
                                                 277, 1234, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3115, 0, 3, 2137, 2173, 949, 2605, 277,
                                                 292, 1279, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3205, 0, 3, 2245, 2305, 1054, 2755, 322,
                                                 343, 1450, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3331, 0, 3, 2305, 2365, 1099, 2845, 343,
                                                 364, 1513, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3457, 0, 3, 2365, 2425, 1144, 2935, 364,
                                                 385, 1576, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3583, 0, 3, 2425, 2485, 1189, 3025, 385,
                                                 406, 1639, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3709, 0, 3, 2485, 2545, 1234, 3115, 406,
                                                 427, 1702, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3835, 3, 469, 472, 1771, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3845, 3, 472, 475, 1777, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3855, 3, 475, 478, 1783, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3865, 3, 478, 481, 1789, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3875, 3, 481, 484, 1795, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3885, 3, 484, 487, 1801, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3895, 3, 487, 490, 1807, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3905, 0, 3, 1765, 3835, 1831, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3935, 0, 3, 1771, 3845, 1849, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3965, 0, 3, 1777, 3855, 1867, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3995, 0, 3, 1783, 3865, 1885, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4025, 0, 3, 1789, 3875, 1903, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4055, 0, 3, 1795, 3885, 1921, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4085, 0, 3, 1801, 3895, 1939, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4115, 0, 3, 1813, 3905, 577, 595, 1993,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4175, 0, 3, 1831, 3935, 595, 613, 2029,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4235, 0, 3, 1849, 3965, 613, 631, 2065,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4295, 0, 3, 1867, 3995, 631, 649, 2101,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4355, 0, 3, 1885, 4025, 649, 667, 2137,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4415, 0, 3, 1903, 4055, 667, 685, 2173,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4475, 0, 3, 1921, 4085, 685, 703, 2209,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4535, 0, 3, 1957, 4115, 739, 769, 2245,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4635, 0, 3, 1993, 4175, 769, 799, 2305,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4735, 0, 3, 2029, 4235, 799, 829, 2365,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4835, 0, 3, 2065, 4295, 829, 859, 2425,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4935, 0, 3, 2101, 4355, 859, 889, 2485,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5035, 0, 3, 2137, 4415, 889, 919, 2545,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5135, 0, 3, 2173, 4475, 919, 949, 2605,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5235, 0, 3, 4115, 4175, 2305, 4735,
                                                 1009, 1054, 2755, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5385, 0, 3, 4175, 4235, 2365, 4835,
                                                 1054, 1099, 2845, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5535, 0, 3, 4235, 4295, 2425, 4935,
                                                 1099, 1144, 2935, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5685, 0, 3, 4295, 4355, 2485, 5035,
                                                 1144, 1189, 3025, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5835, 0, 3, 4355, 4415, 2545, 5135,
                                                 1189, 1234, 3115, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5985, 0, 3, 4535, 4635, 2665, 5235,
                                                 1324, 1387, 3205, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6195, 0, 3, 4635, 4735, 2755, 5385,
                                                 1387, 1450, 3331, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6405, 0, 3, 4735, 4835, 2845, 5535,
                                                 1450, 1513, 3457, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6615, 0, 3, 4835, 4935, 2935, 5685,
                                                 1513, 1576, 3583, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6825, 0, 3, 4935, 5035, 3025, 5835,
                                                 1576, 1639, 3709, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7035, 3, 1765, 1771, 3845, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7050, 3, 1771, 1777, 3855, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7065, 3, 1777, 1783, 3865, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7080, 3, 1783, 1789, 3875, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7095, 3, 1789, 1795, 3885, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7110, 3, 1795, 1801, 3895, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 7125, 0, 3, 3835, 7035, 3935, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7170, 0, 3, 3845, 7050, 3965, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7215, 0, 3, 3855, 7065, 3995, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7260, 0, 3, 3865, 7080, 4025, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7305, 0, 3, 3875, 7095, 4055, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7350, 0, 3, 3885, 7110, 4085, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 7395, 0, 3, 3905, 7125, 1957, 1993,
                                                 4175, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7485, 0, 3, 3935, 7170, 1993, 2029,
                                                 4235, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7575, 0, 3, 3965, 7215, 2029, 2065,
                                                 4295, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7665, 0, 3, 3995, 7260, 2065, 2101,
                                                 4355, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7755, 0, 3, 4025, 7305, 2101, 2137,
                                                 4415, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7845, 0, 3, 4055, 7350, 2137, 2173,
                                                 4475, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7935, 0, 3, 4175, 7485, 2245, 2305,
                                                 4735, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8085, 0, 3, 4235, 7575, 2305, 2365,
                                                 4835, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8235, 0, 3, 4295, 7665, 2365, 2425,
                                                 4935, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8385, 0, 3, 4355, 7755, 2425, 2485,
                                                 5035, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8535, 0, 3, 4415, 7845, 2485, 2545,
                                                 5135, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 8685, 0, 3, 7395, 7485, 4735, 8085,
                                                 2665, 2755, 5385, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 8910, 0, 3, 7485, 7575, 4835, 8235,
                                                 2755, 2845, 5535, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9135, 0, 3, 7575, 7665, 4935, 8385,
                                                 2845, 2935, 5685, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9360, 0, 3, 7665, 7755, 5035, 8535,
                                                 2935, 3025, 5835, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 9585, 0, 3, 7935, 8085, 5385, 8910,
                                                 3205, 3331, 6405, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 9900, 0, 3, 8085, 8235, 5535, 9135,
                                                 3331, 3457, 6615, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 10215, 0, 3, 8235, 8385, 5685, 9360,
                                                 3457, 3583, 6825, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10530, 3, 3835, 3845, 7050, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10551, 3, 3845, 3855, 7065, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10572, 3, 3855, 3865, 7080, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10593, 3, 3865, 3875, 7095, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10614, 3, 3875, 3885, 7110, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 10635, 0, 3, 7035, 10530, 7170, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10698, 0, 3, 7050, 10551, 7215, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10761, 0, 3, 7065, 10572, 7260, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10824, 0, 3, 7080, 10593, 7305, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10887, 0, 3, 7095, 10614, 7350, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 10950, 0, 3, 7125, 10635, 4115, 4175,
                                                 7485, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11076, 0, 3, 7170, 10698, 4175, 4235,
                                                 7575, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11202, 0, 3, 7215, 10761, 4235, 4295,
                                                 7665, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11328, 0, 3, 7260, 10824, 4295, 4355,
                                                 7755, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11454, 0, 3, 7305, 10887, 4355, 4415,
                                                 7845, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 11580, 0, 3, 7395, 10950, 4535, 4635,
                                                 7935, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 11790, 0, 3, 7485, 11076, 4635, 4735,
                                                 8085, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12000, 0, 3, 7575, 11202, 4735, 4835,
                                                 8235, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12210, 0, 3, 7665, 11328, 4835, 4935,
                                                 8385, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12420, 0, 3, 7755, 11454, 4935, 5035,
                                                 8535, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 12630, 0, 3, 10950, 11076, 8085, 12000,
                                                 5235, 5385, 8910, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 12945, 0, 3, 11076, 11202, 8235, 12210,
                                                 5385, 5535, 9135, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 13260, 0, 3, 11202, 11328, 8385, 12420,
                                                 5535, 5685, 9360, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 13575, 0, 3, 11580, 11790, 8685, 12630,
                                                 5985, 6195, 9585, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 14016, 0, 3, 11790, 12000, 8910, 12945,
                                                 6195, 6405, 9900, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 14457, 0, 3, 12000, 12210, 9135, 13260,
                                                 6405, 6615, 10215, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 14898, 3, 7035, 7050, 10551, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 14926, 3, 7050, 7065, 10572, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 14954, 3, 7065, 7080, 10593, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 14982, 3, 7080, 7095, 10614, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 15010, 0, 3, 10530, 14898, 10698, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 15094, 0, 3, 10551, 14926, 10761, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 15178, 0, 3, 10572, 14954, 10824, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 15262, 0, 3, 10593, 14982, 10887, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 15346, 0, 3, 10635, 15010, 7395, 7485,
                                                 11076, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 15514, 0, 3, 10698, 15094, 7485, 7575,
                                                 11202, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 15682, 0, 3, 10761, 15178, 7575, 7665,
                                                 11328, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 15850, 0, 3, 10824, 15262, 7665, 7755,
                                                 11454, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 16018, 0, 3, 11076, 15514, 7935, 8085,
                                                 12000, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 16298, 0, 3, 11202, 15682, 8085, 8235,
                                                 12210, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 16578, 0, 3, 11328, 15850, 8235, 8385,
                                                 12420, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 16858, 0, 3, 15346, 15514, 12000, 16298,
                                                 8685, 8910, 12945, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 17278, 0, 3, 15514, 15682, 12210, 16578,
                                                 8910, 9135, 13260, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 17698, 0, 3, 16018, 16298, 12945, 17278,
                                                 9585, 9900, 14457, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 18286, 3, 10530, 10551, 14926, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 18322, 3, 10551, 10572, 14954, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 18358, 3, 10572, 10593, 14982, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 18394, 0, 3, 14898, 18286, 15094, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 18502, 0, 3, 14926, 18322, 15178, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 18610, 0, 3, 14954, 18358, 15262, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 18718, 0, 3, 15010, 18394, 10950, 11076,
                                                 15514, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 18934, 0, 3, 15094, 18502, 11076, 11202,
                                                 15682, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 19150, 0, 3, 15178, 18610, 11202, 11328,
                                                 15850, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 19366, 0, 3, 15346, 18718, 11580, 11790,
                                                 16018, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 19726, 0, 3, 15514, 18934, 11790, 12000,
                                                 16298, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 20086, 0, 3, 15682, 19150, 12000, 12210,
                                                 16578, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 20446, 0, 3, 18718, 18934, 16298, 20086,
                                                 12630, 12945, 17278, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 20986, 0, 3, 19366, 19726, 16858, 20446,
                                                 13575, 14016, 17698, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 21742, 20986, 756, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 22498, 21742, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 22498, 15, nmax);
}

}  // namespace simdt2ceri
