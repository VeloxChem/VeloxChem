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


#include "SimdElectronRepulsionGeom10RsRecGH.hpp"

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
#include "SimdGeometryG1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_gh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_gh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 21191, 19136, 1890, nvalues);

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

            compute_prim_gs_electron_repulsion_0(buffer, 318, 0, 82, 88, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 333, 0, 88, 94, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 348, 0, 94, 100, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 363, 0, 100, 106, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 378, 0, 106, 112, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 393, 0, 112, 118, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 408, 0, 130, 136, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 423, 0, 136, 142, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 438, 0, 142, 148, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 453, 0, 148, 154, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 468, 0, 154, 160, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 483, 0, 160, 166, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 498, 0, 178, 188, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 519, 0, 188, 198, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 540, 0, 198, 208, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 561, 0, 208, 218, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 582, 0, 218, 228, 393, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 603, 0, 248, 258, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 624, 0, 258, 268, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 645, 0, 268, 278, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 666, 0, 278, 288, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 687, 0, 288, 298, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 708, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 711, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 714, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 717, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 720, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 723, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 726, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 729, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 732, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 735, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 738, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 741, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 744, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 747, 3, 27, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 750, 3, 9, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 759, 3, 10, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 768, 3, 11, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 777, 3, 12, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 786, 3, 13, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 795, 3, 14, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 804, 3, 15, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 813, 3, 20, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 822, 3, 21, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 831, 3, 22, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 840, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 849, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 858, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 867, 3, 26, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 876, 0, 3, 31, 750, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 894, 0, 3, 34, 759, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 912, 0, 3, 37, 768, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 930, 0, 3, 40, 777, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 948, 0, 3, 43, 786, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 966, 0, 3, 46, 795, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 984, 0, 3, 49, 804, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1002, 0, 3, 58, 813, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1020, 0, 3, 61, 822, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1038, 0, 3, 64, 831, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1056, 0, 3, 67, 840, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1074, 0, 3, 70, 849, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1092, 0, 3, 73, 858, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1110, 0, 3, 76, 867, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1128, 0, 3, 82, 876, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1158, 0, 3, 88, 894, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1188, 0, 3, 94, 912, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1218, 0, 3, 100, 930, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1248, 0, 3, 106, 948, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1278, 0, 3, 112, 966, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1308, 0, 3, 118, 984, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1338, 0, 3, 130, 1002, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1368, 0, 3, 136, 1020, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1398, 0, 3, 142, 1038, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1428, 0, 3, 148, 1056, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1458, 0, 3, 154, 1074, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1488, 0, 3, 160, 1092, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1518, 0, 3, 166, 1110, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1548, 0, 3, 188, 1188, 333, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1593, 0, 3, 198, 1218, 348, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1638, 0, 3, 208, 1248, 363, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1683, 0, 3, 218, 1278, 378, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1728, 0, 3, 228, 1308, 393, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1773, 0, 3, 258, 1398, 423, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1818, 0, 3, 268, 1428, 438, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1863, 0, 3, 278, 1458, 453, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1908, 0, 3, 288, 1488, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1953, 0, 3, 298, 1518, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1998, 0, 3, 318, 1548, 498, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2061, 0, 3, 333, 1593, 519, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2124, 0, 3, 348, 1638, 540, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2187, 0, 3, 363, 1683, 561, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2250, 0, 3, 378, 1728, 582, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2313, 0, 3, 408, 1773, 603, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2376, 0, 3, 423, 1818, 624, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2439, 0, 3, 438, 1863, 645, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2502, 0, 3, 453, 1908, 666, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2565, 0, 3, 468, 1953, 687, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2628, 3, 9, 10, 711, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2634, 3, 10, 11, 714, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2640, 3, 11, 12, 717, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2646, 3, 12, 13, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2652, 3, 13, 14, 723, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2658, 3, 14, 15, 726, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2664, 3, 20, 21, 732, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2670, 3, 21, 22, 735, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2676, 3, 22, 23, 738, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2682, 3, 23, 24, 741, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2688, 3, 24, 25, 744, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2694, 3, 25, 26, 747, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2700, 0, 3, 708, 2628, 759, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2718, 0, 3, 711, 2634, 768, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2736, 0, 3, 714, 2640, 777, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2754, 0, 3, 717, 2646, 786, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2772, 0, 3, 720, 2652, 795, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2790, 0, 3, 723, 2658, 804, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2808, 0, 3, 729, 2664, 822, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2826, 0, 3, 732, 2670, 831, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2844, 0, 3, 735, 2676, 840, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2862, 0, 3, 738, 2682, 849, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2880, 0, 3, 741, 2688, 858, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2898, 0, 3, 744, 2694, 867, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2916, 0, 3, 750, 2700, 82, 88, 894,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2952, 0, 3, 759, 2718, 88, 94, 912,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2988, 0, 3, 768, 2736, 94, 100, 930,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3024, 0, 3, 777, 2754, 100, 106, 948,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3060, 0, 3, 786, 2772, 106, 112, 966,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3096, 0, 3, 795, 2790, 112, 118, 984,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3132, 0, 3, 813, 2808, 130, 136, 1020,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3168, 0, 3, 822, 2826, 136, 142, 1038,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3204, 0, 3, 831, 2844, 142, 148, 1056,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3240, 0, 3, 840, 2862, 148, 154, 1074,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3276, 0, 3, 849, 2880, 154, 160, 1092,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3312, 0, 3, 858, 2898, 160, 166, 1110,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3348, 0, 3, 894, 2952, 178, 188, 1188,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3408, 0, 3, 912, 2988, 188, 198, 1218,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3468, 0, 3, 930, 3024, 198, 208, 1248,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3528, 0, 3, 948, 3060, 208, 218, 1278,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3588, 0, 3, 966, 3096, 218, 228, 1308,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3648, 0, 3, 1020, 3168, 248, 258, 1398,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3708, 0, 3, 1038, 3204, 258, 268, 1428,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3768, 0, 3, 1056, 3240, 268, 278, 1458,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3828, 0, 3, 1074, 3276, 278, 288, 1488,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3888, 0, 3, 1092, 3312, 288, 298, 1518,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3948, 0, 3, 2916, 2952, 1188, 3408, 318,
                                                 333, 1593, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4038, 0, 3, 2952, 2988, 1218, 3468, 333,
                                                 348, 1638, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4128, 0, 3, 2988, 3024, 1248, 3528, 348,
                                                 363, 1683, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4218, 0, 3, 3024, 3060, 1278, 3588, 363,
                                                 378, 1728, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4308, 0, 3, 3132, 3168, 1398, 3708, 408,
                                                 423, 1818, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4398, 0, 3, 3168, 3204, 1428, 3768, 423,
                                                 438, 1863, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4488, 0, 3, 3204, 3240, 1458, 3828, 438,
                                                 453, 1908, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4578, 0, 3, 3240, 3276, 1488, 3888, 453,
                                                 468, 1953, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4668, 0, 3, 3348, 3408, 1593, 4038, 498,
                                                 519, 2124, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4794, 0, 3, 3408, 3468, 1638, 4128, 519,
                                                 540, 2187, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4920, 0, 3, 3468, 3528, 1683, 4218, 540,
                                                 561, 2250, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5046, 0, 3, 3648, 3708, 1818, 4398, 603,
                                                 624, 2439, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5172, 0, 3, 3708, 3768, 1863, 4488, 624,
                                                 645, 2502, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5298, 0, 3, 3768, 3828, 1908, 4578, 645,
                                                 666, 2565, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5424, 3, 708, 711, 2634, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5434, 3, 711, 714, 2640, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5444, 3, 714, 717, 2646, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5454, 3, 717, 720, 2652, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5464, 3, 720, 723, 2658, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5474, 3, 729, 732, 2670, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5484, 3, 732, 735, 2676, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5494, 3, 735, 738, 2682, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5504, 3, 738, 741, 2688, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5514, 3, 741, 744, 2694, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 5524, 0, 3, 2628, 5424, 2718, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5554, 0, 3, 2634, 5434, 2736, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5584, 0, 3, 2640, 5444, 2754, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5614, 0, 3, 2646, 5454, 2772, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5644, 0, 3, 2652, 5464, 2790, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5674, 0, 3, 2664, 5474, 2826, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5704, 0, 3, 2670, 5484, 2844, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5734, 0, 3, 2676, 5494, 2862, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5764, 0, 3, 2682, 5504, 2880, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5794, 0, 3, 2688, 5514, 2898, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 5824, 0, 3, 2700, 5524, 876, 894, 2952,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5884, 0, 3, 2718, 5554, 894, 912, 2988,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5944, 0, 3, 2736, 5584, 912, 930, 3024,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6004, 0, 3, 2754, 5614, 930, 948, 3060,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6064, 0, 3, 2772, 5644, 948, 966, 3096,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6124, 0, 3, 2808, 5674, 1002, 1020,
                                                 3168, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6184, 0, 3, 2826, 5704, 1020, 1038,
                                                 3204, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6244, 0, 3, 2844, 5734, 1038, 1056,
                                                 3240, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6304, 0, 3, 2862, 5764, 1056, 1074,
                                                 3276, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6364, 0, 3, 2880, 5794, 1074, 1092,
                                                 3312, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6424, 0, 3, 2916, 5824, 1128, 1158,
                                                 3348, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6524, 0, 3, 2952, 5884, 1158, 1188,
                                                 3408, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6624, 0, 3, 2988, 5944, 1188, 1218,
                                                 3468, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6724, 0, 3, 3024, 6004, 1218, 1248,
                                                 3528, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6824, 0, 3, 3060, 6064, 1248, 1278,
                                                 3588, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6924, 0, 3, 3132, 6124, 1338, 1368,
                                                 3648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7024, 0, 3, 3168, 6184, 1368, 1398,
                                                 3708, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7124, 0, 3, 3204, 6244, 1398, 1428,
                                                 3768, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7224, 0, 3, 3240, 6304, 1428, 1458,
                                                 3828, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7324, 0, 3, 3276, 6364, 1458, 1488,
                                                 3888, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7424, 0, 3, 5824, 5884, 3408, 6624,
                                                 1548, 1593, 4038, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7574, 0, 3, 5884, 5944, 3468, 6724,
                                                 1593, 1638, 4128, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7724, 0, 3, 5944, 6004, 3528, 6824,
                                                 1638, 1683, 4218, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7874, 0, 3, 6124, 6184, 3708, 7124,
                                                 1773, 1818, 4398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8024, 0, 3, 6184, 6244, 3768, 7224,
                                                 1818, 1863, 4488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8174, 0, 3, 6244, 6304, 3828, 7324,
                                                 1863, 1908, 4578, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8324, 0, 3, 6424, 6524, 3948, 7424,
                                                 1998, 2061, 4668, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8534, 0, 3, 6524, 6624, 4038, 7574,
                                                 2061, 2124, 4794, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8744, 0, 3, 6624, 6724, 4128, 7724,
                                                 2124, 2187, 4920, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8954, 0, 3, 6924, 7024, 4308, 7874,
                                                 2313, 2376, 5046, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9164, 0, 3, 7024, 7124, 4398, 8024,
                                                 2376, 2439, 5172, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9374, 0, 3, 7124, 7224, 4488, 8174,
                                                 2439, 2502, 5298, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9584, 3, 2628, 2634, 5434, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9599, 3, 2634, 2640, 5444, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9614, 3, 2640, 2646, 5454, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9629, 3, 2646, 2652, 5464, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9644, 3, 2664, 2670, 5484, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9659, 3, 2670, 2676, 5494, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9674, 3, 2676, 2682, 5504, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9689, 3, 2682, 2688, 5514, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 9704, 0, 3, 5424, 9584, 5554, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9749, 0, 3, 5434, 9599, 5584, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9794, 0, 3, 5444, 9614, 5614, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9839, 0, 3, 5454, 9629, 5644, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9884, 0, 3, 5474, 9644, 5704, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9929, 0, 3, 5484, 9659, 5734, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9974, 0, 3, 5494, 9674, 5764, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10019, 0, 3, 5504, 9689, 5794, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 10064, 0, 3, 5524, 9704, 2916, 2952,
                                                 5884, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10154, 0, 3, 5554, 9749, 2952, 2988,
                                                 5944, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10244, 0, 3, 5584, 9794, 2988, 3024,
                                                 6004, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10334, 0, 3, 5614, 9839, 3024, 3060,
                                                 6064, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10424, 0, 3, 5674, 9884, 3132, 3168,
                                                 6184, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10514, 0, 3, 5704, 9929, 3168, 3204,
                                                 6244, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10604, 0, 3, 5734, 9974, 3204, 3240,
                                                 6304, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10694, 0, 3, 5764, 10019, 3240, 3276,
                                                 6364, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10784, 0, 3, 5884, 10154, 3348, 3408,
                                                 6624, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 10934, 0, 3, 5944, 10244, 3408, 3468,
                                                 6724, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11084, 0, 3, 6004, 10334, 3468, 3528,
                                                 6824, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11234, 0, 3, 6184, 10514, 3648, 3708,
                                                 7124, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11384, 0, 3, 6244, 10604, 3708, 3768,
                                                 7224, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11534, 0, 3, 6304, 10694, 3768, 3828,
                                                 7324, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11684, 0, 3, 10064, 10154, 6624, 10934,
                                                 3948, 4038, 7574, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 11909, 0, 3, 10154, 10244, 6724, 11084,
                                                 4038, 4128, 7724, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 12134, 0, 3, 10424, 10514, 7124, 11384,
                                                 4308, 4398, 8024, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 12359, 0, 3, 10514, 10604, 7224, 11534,
                                                 4398, 4488, 8174, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12584, 0, 3, 10784, 10934, 7574, 11909,
                                                 4668, 4794, 8744, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 12899, 0, 3, 11234, 11384, 8024, 12359,
                                                 5046, 5172, 9374, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13214, 3, 5424, 5434, 9599, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13235, 3, 5434, 5444, 9614, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13256, 3, 5444, 5454, 9629, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13277, 3, 5474, 5484, 9659, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13298, 3, 5484, 5494, 9674, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 13319, 3, 5494, 5504, 9689, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 13340, 0, 3, 9584, 13214, 9749, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13403, 0, 3, 9599, 13235, 9794, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13466, 0, 3, 9614, 13256, 9839, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13529, 0, 3, 9644, 13277, 9929, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13592, 0, 3, 9659, 13298, 9974, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 13655, 0, 3, 9674, 13319, 10019, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 13718, 0, 3, 9704, 13340, 5824, 5884,
                                                 10154, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13844, 0, 3, 9749, 13403, 5884, 5944,
                                                 10244, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 13970, 0, 3, 9794, 13466, 5944, 6004,
                                                 10334, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14096, 0, 3, 9884, 13529, 6124, 6184,
                                                 10514, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14222, 0, 3, 9929, 13592, 6184, 6244,
                                                 10604, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 14348, 0, 3, 9974, 13655, 6244, 6304,
                                                 10694, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14474, 0, 3, 10064, 13718, 6424, 6524,
                                                 10784, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14684, 0, 3, 10154, 13844, 6524, 6624,
                                                 10934, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 14894, 0, 3, 10244, 13970, 6624, 6724,
                                                 11084, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15104, 0, 3, 10424, 14096, 6924, 7024,
                                                 11234, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15314, 0, 3, 10514, 14222, 7024, 7124,
                                                 11384, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 15524, 0, 3, 10604, 14348, 7124, 7224,
                                                 11534, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 15734, 0, 3, 13718, 13844, 10934, 14894,
                                                 7424, 7574, 11909, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 16049, 0, 3, 14096, 14222, 11384, 15524,
                                                 7874, 8024, 12359, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 16364, 0, 3, 14474, 14684, 11684, 15734,
                                                 8324, 8534, 12584, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 16805, 0, 3, 15104, 15314, 12134, 16049,
                                                 8954, 9164, 12899, ncols, alpha, beta, p);

            simdgeo::geom_g_x(buffer, 17246, 15104, 16805, 1, 21, ncols, alpha);

            simdgeo::geom_g_y(buffer, 17561, 15104, 16805, 1, 21, ncols, alpha);

            simdgeo::geom_g_z(buffer, 17876, 15104, 16805, 1, 21, ncols, alpha);

            simdgeo::geom_g_x(buffer, 18191, 14474, 16364, 1, 21, ncols, alpha);

            simdgeo::geom_g_y(buffer, 18506, 14474, 16364, 1, 21, ncols, alpha);

            simdgeo::geom_g_z(buffer, 18821, 14474, 16364, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 19136, 18191, 945, ncols);

            simdfunc::contract_primitives(buffer, 20081, 17246, 945, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 21026, 20081, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 21026, 11, nmax);

    simdtrf::transform_h_inner(buffer, 21026, 20396, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 99 * nvalues, nvalues, buffer, 21026, 11, nmax);

    simdtrf::transform_h_inner(buffer, 21026, 20711, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 198 * nvalues, nvalues, buffer, 21026, 11, nmax);

    simdtrf::transform_h_inner(buffer, 21026, 19136, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 297 * nvalues, nvalues, buffer, 21026, 11, nmax);

    simdtrf::transform_h_inner(buffer, 21026, 19451, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 396 * nvalues, nvalues, buffer, 21026, 11, nmax);

    simdtrf::transform_h_inner(buffer, 21026, 19766, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 495 * nvalues, nvalues, buffer, 21026, 11, nmax);
}

}  // namespace simdt2ceri
