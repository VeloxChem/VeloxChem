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


#include "SimdElectronRepulsionRsRecHG.hpp"

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
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_hg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 10677, 9858, 630, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 7, 8, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 8, 9, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 9, 10, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 10, 11, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 11, 12, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 12, 13, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 13, 14, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 17, 18, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 18, 19, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 19, 20, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 20, 21, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 21, 22, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 22, 23, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 44, 47, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 53, 56, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 56, 59, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 59, 62, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 62, 65, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 65, 68, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 68, 71, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 71, 74, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 304, 0, 80, 86, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 319, 0, 86, 92, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 334, 0, 92, 98, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 349, 0, 98, 104, 214, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 364, 0, 104, 110, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 379, 0, 122, 128, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 394, 0, 128, 134, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 409, 0, 134, 140, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 424, 0, 140, 146, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 439, 0, 146, 152, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 454, 0, 164, 174, 304, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 475, 0, 174, 184, 319, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 496, 0, 184, 194, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 517, 0, 194, 204, 349, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 538, 0, 204, 214, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 559, 0, 234, 244, 379, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 580, 0, 244, 254, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 601, 0, 254, 264, 409, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 622, 0, 264, 274, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 643, 0, 274, 284, 439, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 664, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 667, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 670, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 673, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 676, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 679, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 682, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 685, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 688, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 691, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 694, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 697, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 700, 3, 9, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 709, 3, 10, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 718, 3, 11, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 727, 3, 12, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 736, 3, 13, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 745, 3, 14, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 754, 3, 19, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 763, 3, 20, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 772, 3, 21, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 781, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 790, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 799, 3, 24, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 808, 0, 3, 32, 700, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 826, 0, 3, 35, 709, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 844, 0, 3, 38, 718, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 862, 0, 3, 41, 727, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 880, 0, 3, 44, 736, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 898, 0, 3, 47, 745, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 916, 0, 3, 59, 754, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 934, 0, 3, 62, 763, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 952, 0, 3, 65, 772, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 970, 0, 3, 68, 781, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 988, 0, 3, 71, 790, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1006, 0, 3, 74, 799, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1024, 0, 3, 86, 826, 184, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1054, 0, 3, 92, 844, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1084, 0, 3, 98, 862, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1114, 0, 3, 104, 880, 214, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1144, 0, 3, 110, 898, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1174, 0, 3, 128, 934, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1204, 0, 3, 134, 952, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1234, 0, 3, 140, 970, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1264, 0, 3, 146, 988, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1294, 0, 3, 152, 1006, 294, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1324, 0, 3, 184, 1054, 319, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1369, 0, 3, 194, 1084, 334, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1414, 0, 3, 204, 1114, 349, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1459, 0, 3, 214, 1144, 364, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1504, 0, 3, 254, 1204, 394, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1549, 0, 3, 264, 1234, 409, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1594, 0, 3, 274, 1264, 424, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1639, 0, 3, 284, 1294, 439, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1684, 0, 3, 319, 1369, 496, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1747, 0, 3, 334, 1414, 517, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1810, 0, 3, 349, 1459, 538, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1873, 0, 3, 394, 1549, 601, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1936, 0, 3, 409, 1594, 622, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1999, 0, 3, 424, 1639, 643, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2062, 3, 9, 10, 667, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2068, 3, 10, 11, 670, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2074, 3, 11, 12, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2080, 3, 12, 13, 676, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2086, 3, 13, 14, 679, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2092, 3, 19, 20, 685, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2098, 3, 20, 21, 688, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2104, 3, 21, 22, 691, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2110, 3, 22, 23, 694, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2116, 3, 23, 24, 697, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2122, 0, 3, 664, 2062, 709, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2140, 0, 3, 667, 2068, 718, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2158, 0, 3, 670, 2074, 727, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2176, 0, 3, 673, 2080, 736, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2194, 0, 3, 676, 2086, 745, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2212, 0, 3, 682, 2092, 763, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2230, 0, 3, 685, 2098, 772, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2248, 0, 3, 688, 2104, 781, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2266, 0, 3, 691, 2110, 790, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2284, 0, 3, 694, 2116, 799, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2302, 0, 3, 700, 2122, 80, 86, 826,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2338, 0, 3, 709, 2140, 86, 92, 844,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2374, 0, 3, 718, 2158, 92, 98, 862,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2410, 0, 3, 727, 2176, 98, 104, 880,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2446, 0, 3, 736, 2194, 104, 110, 898,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2482, 0, 3, 754, 2212, 122, 128, 934,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2518, 0, 3, 763, 2230, 128, 134, 952,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2554, 0, 3, 772, 2248, 134, 140, 970,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2590, 0, 3, 781, 2266, 140, 146, 988,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2626, 0, 3, 790, 2284, 146, 152, 1006,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2662, 0, 3, 808, 2302, 164, 174, 1024,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2722, 0, 3, 826, 2338, 174, 184, 1054,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2782, 0, 3, 844, 2374, 184, 194, 1084,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2842, 0, 3, 862, 2410, 194, 204, 1114,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2902, 0, 3, 880, 2446, 204, 214, 1144,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2962, 0, 3, 916, 2482, 234, 244, 1174,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3022, 0, 3, 934, 2518, 244, 254, 1204,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3082, 0, 3, 952, 2554, 254, 264, 1234,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3142, 0, 3, 970, 2590, 264, 274, 1264,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3202, 0, 3, 988, 2626, 274, 284, 1294,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3262, 0, 3, 2302, 2338, 1054, 2782, 304,
                                                 319, 1369, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3352, 0, 3, 2338, 2374, 1084, 2842, 319,
                                                 334, 1414, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3442, 0, 3, 2374, 2410, 1114, 2902, 334,
                                                 349, 1459, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3532, 0, 3, 2482, 2518, 1204, 3082, 379,
                                                 394, 1549, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3622, 0, 3, 2518, 2554, 1234, 3142, 394,
                                                 409, 1594, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3712, 0, 3, 2554, 2590, 1264, 3202, 409,
                                                 424, 1639, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3802, 0, 3, 2662, 2722, 1324, 3262, 454,
                                                 475, 1684, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3928, 0, 3, 2722, 2782, 1369, 3352, 475,
                                                 496, 1747, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4054, 0, 3, 2782, 2842, 1414, 3442, 496,
                                                 517, 1810, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4180, 0, 3, 2962, 3022, 1504, 3532, 559,
                                                 580, 1873, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4306, 0, 3, 3022, 3082, 1549, 3622, 580,
                                                 601, 1936, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4432, 0, 3, 3082, 3142, 1594, 3712, 601,
                                                 622, 1999, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4558, 3, 664, 667, 2068, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4568, 3, 667, 670, 2074, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4578, 3, 670, 673, 2080, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4588, 3, 673, 676, 2086, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4598, 3, 682, 685, 2098, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4608, 3, 685, 688, 2104, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4618, 3, 688, 691, 2110, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4628, 3, 691, 694, 2116, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4638, 0, 3, 2062, 4558, 2140, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4668, 0, 3, 2068, 4568, 2158, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4698, 0, 3, 2074, 4578, 2176, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4728, 0, 3, 2080, 4588, 2194, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4758, 0, 3, 2092, 4598, 2230, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4788, 0, 3, 2098, 4608, 2248, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4818, 0, 3, 2104, 4618, 2266, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4848, 0, 3, 2110, 4628, 2284, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4878, 0, 3, 2122, 4638, 808, 826, 2338,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4938, 0, 3, 2140, 4668, 826, 844, 2374,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4998, 0, 3, 2158, 4698, 844, 862, 2410,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5058, 0, 3, 2176, 4728, 862, 880, 2446,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5118, 0, 3, 2212, 4758, 916, 934, 2518,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5178, 0, 3, 2230, 4788, 934, 952, 2554,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5238, 0, 3, 2248, 4818, 952, 970, 2590,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5298, 0, 3, 2266, 4848, 970, 988, 2626,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5358, 0, 3, 2338, 4938, 1024, 1054,
                                                 2782, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5458, 0, 3, 2374, 4998, 1054, 1084,
                                                 2842, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5558, 0, 3, 2410, 5058, 1084, 1114,
                                                 2902, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5658, 0, 3, 2518, 5178, 1174, 1204,
                                                 3082, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5758, 0, 3, 2554, 5238, 1204, 1234,
                                                 3142, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5858, 0, 3, 2590, 5298, 1234, 1264,
                                                 3202, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 5958, 0, 3, 4878, 4938, 2782, 5458,
                                                 1324, 1369, 3352, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6108, 0, 3, 4938, 4998, 2842, 5558,
                                                 1369, 1414, 3442, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6258, 0, 3, 5118, 5178, 3082, 5758,
                                                 1504, 1549, 3622, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6408, 0, 3, 5178, 5238, 3142, 5858,
                                                 1549, 1594, 3712, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6558, 0, 3, 5358, 5458, 3352, 6108,
                                                 1684, 1747, 4054, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 6768, 0, 3, 5658, 5758, 3622, 6408,
                                                 1873, 1936, 4432, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6978, 3, 2062, 2068, 4568, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6993, 3, 2068, 2074, 4578, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7008, 3, 2074, 2080, 4588, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7023, 3, 2092, 2098, 4608, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7038, 3, 2098, 2104, 4618, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7053, 3, 2104, 2110, 4628, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 7068, 0, 3, 4558, 6978, 4668, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7113, 0, 3, 4568, 6993, 4698, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7158, 0, 3, 4578, 7008, 4728, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7203, 0, 3, 4598, 7023, 4788, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7248, 0, 3, 4608, 7038, 4818, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7293, 0, 3, 4618, 7053, 4848, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 7338, 0, 3, 4638, 7068, 2302, 2338,
                                                 4938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7428, 0, 3, 4668, 7113, 2338, 2374,
                                                 4998, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7518, 0, 3, 4698, 7158, 2374, 2410,
                                                 5058, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7608, 0, 3, 4758, 7203, 2482, 2518,
                                                 5178, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7698, 0, 3, 4788, 7248, 2518, 2554,
                                                 5238, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7788, 0, 3, 4818, 7293, 2554, 2590,
                                                 5298, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7878, 0, 3, 4878, 7338, 2662, 2722,
                                                 5358, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8028, 0, 3, 4938, 7428, 2722, 2782,
                                                 5458, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8178, 0, 3, 4998, 7518, 2782, 2842,
                                                 5558, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8328, 0, 3, 5118, 7608, 2962, 3022,
                                                 5658, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8478, 0, 3, 5178, 7698, 3022, 3082,
                                                 5758, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8628, 0, 3, 5238, 7788, 3082, 3142,
                                                 5858, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 8778, 0, 3, 7338, 7428, 5458, 8178,
                                                 3262, 3352, 6108, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9003, 0, 3, 7608, 7698, 5758, 8628,
                                                 3532, 3622, 6408, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 9228, 0, 3, 7878, 8028, 5958, 8778,
                                                 3802, 3928, 6558, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 9543, 0, 3, 8328, 8478, 6258, 9003,
                                                 4180, 4306, 6768, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9858, 9228, 630, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 10488, 10173, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 10488, 9, nmax);

    simdtrf::transform_g_inner(buffer, 10488, 9858, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 99 * nvalues, nvalues, buffer, 10488, 9, nmax);
}

}  // namespace simdt2ceri
