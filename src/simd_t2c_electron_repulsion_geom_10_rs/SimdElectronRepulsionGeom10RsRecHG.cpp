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


#include "SimdElectronRepulsionGeom10RsRecHG.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_hg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_hg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 20767, 18688, 1890, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 18, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 7, 8, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 8, 9, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 9, 10, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 10, 11, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 11, 12, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 12, 13, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 13, 14, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 14, 15, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 15, 16, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 22, 23, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 23, 24, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 24, 25, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 25, 26, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 26, 27, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 27, 28, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 262, 0, 57, 60, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 60, 63, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 282, 0, 63, 66, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 292, 0, 66, 69, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 302, 0, 69, 72, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 72, 75, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 75, 78, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 332, 0, 84, 90, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 347, 0, 90, 96, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 362, 0, 96, 102, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 377, 0, 102, 108, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 392, 0, 108, 114, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 407, 0, 114, 120, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 422, 0, 120, 126, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 437, 0, 138, 144, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 452, 0, 144, 150, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 467, 0, 150, 156, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 482, 0, 156, 162, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 497, 0, 162, 168, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 512, 0, 168, 174, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 527, 0, 174, 180, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 542, 0, 192, 202, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 563, 0, 202, 212, 377, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 584, 0, 212, 222, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 605, 0, 222, 232, 407, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 626, 0, 232, 242, 422, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 647, 0, 262, 272, 467, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 668, 0, 272, 282, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 689, 0, 282, 292, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 710, 0, 292, 302, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 731, 0, 302, 312, 527, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 752, 0, 332, 347, 542, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 780, 0, 347, 362, 563, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 808, 0, 362, 377, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 836, 0, 377, 392, 605, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 864, 0, 392, 407, 626, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 892, 0, 437, 452, 647, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 920, 0, 452, 467, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 948, 0, 467, 482, 689, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 976, 0, 482, 497, 710, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1004, 0, 497, 512, 731, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1032, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1035, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1038, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1041, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1044, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1047, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1050, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1053, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1056, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1059, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1062, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1065, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1068, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1071, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1074, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1077, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1080, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1089, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1098, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1107, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1116, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1125, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1134, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1143, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1152, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1161, 3, 22, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1170, 3, 23, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1179, 3, 24, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1188, 3, 25, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1197, 3, 26, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1206, 3, 27, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1215, 3, 28, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1224, 0, 3, 33, 1089, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1242, 0, 3, 36, 1098, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1260, 0, 3, 39, 1107, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1278, 0, 3, 42, 1116, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1296, 0, 3, 45, 1125, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1314, 0, 3, 48, 1134, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1332, 0, 3, 51, 1143, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1350, 0, 3, 60, 1161, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1368, 0, 3, 63, 1170, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1386, 0, 3, 66, 1179, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1404, 0, 3, 69, 1188, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1422, 0, 3, 72, 1197, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1440, 0, 3, 75, 1206, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1458, 0, 3, 78, 1215, 186, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1476, 0, 3, 96, 1242, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1506, 0, 3, 102, 1260, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1536, 0, 3, 108, 1278, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1566, 0, 3, 114, 1296, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1596, 0, 3, 120, 1314, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1626, 0, 3, 126, 1332, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1656, 0, 3, 150, 1368, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1686, 0, 3, 156, 1386, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1716, 0, 3, 162, 1404, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1746, 0, 3, 168, 1422, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1776, 0, 3, 174, 1440, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1806, 0, 3, 180, 1458, 322, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1836, 0, 3, 202, 1506, 362, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1881, 0, 3, 212, 1536, 377, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1926, 0, 3, 222, 1566, 392, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1971, 0, 3, 232, 1596, 407, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2016, 0, 3, 242, 1626, 422, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2061, 0, 3, 272, 1686, 467, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2106, 0, 3, 282, 1716, 482, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2151, 0, 3, 292, 1746, 497, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2196, 0, 3, 302, 1776, 512, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2241, 0, 3, 312, 1806, 527, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2286, 0, 3, 362, 1881, 563, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2349, 0, 3, 377, 1926, 584, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2412, 0, 3, 392, 1971, 605, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2475, 0, 3, 407, 2016, 626, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2538, 0, 3, 467, 2106, 668, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2601, 0, 3, 482, 2151, 689, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2664, 0, 3, 497, 2196, 710, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2727, 0, 3, 512, 2241, 731, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2790, 0, 3, 563, 2349, 808, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2874, 0, 3, 584, 2412, 836, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2958, 0, 3, 605, 2475, 864, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3042, 0, 3, 668, 2601, 948, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3126, 0, 3, 689, 2664, 976, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3210, 0, 3, 710, 2727, 1004, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3294, 3, 9, 10, 1035, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3300, 3, 10, 11, 1038, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3306, 3, 11, 12, 1041, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3312, 3, 12, 13, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3318, 3, 13, 14, 1047, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3324, 3, 14, 15, 1050, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3330, 3, 15, 16, 1053, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3336, 3, 21, 22, 1059, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3342, 3, 22, 23, 1062, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3348, 3, 23, 24, 1065, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3354, 3, 24, 25, 1068, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3360, 3, 25, 26, 1071, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3366, 3, 26, 27, 1074, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3372, 3, 27, 28, 1077, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3378, 0, 3, 1032, 3294, 1089, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3396, 0, 3, 1035, 3300, 1098, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3414, 0, 3, 1038, 3306, 1107, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3432, 0, 3, 1041, 3312, 1116, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3450, 0, 3, 1044, 3318, 1125, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3468, 0, 3, 1047, 3324, 1134, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3486, 0, 3, 1050, 3330, 1143, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3504, 0, 3, 1056, 3336, 1161, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3522, 0, 3, 1059, 3342, 1170, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3540, 0, 3, 1062, 3348, 1179, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3558, 0, 3, 1065, 3354, 1188, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3576, 0, 3, 1068, 3360, 1197, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3594, 0, 3, 1071, 3366, 1206, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3612, 0, 3, 1074, 3372, 1215, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 3630, 0, 3, 1080, 3378, 84, 90, 1224,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3666, 0, 3, 1089, 3396, 90, 96, 1242,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3702, 0, 3, 1098, 3414, 96, 102, 1260,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3738, 0, 3, 1107, 3432, 102, 108, 1278,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3774, 0, 3, 1116, 3450, 108, 114, 1296,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3810, 0, 3, 1125, 3468, 114, 120, 1314,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3846, 0, 3, 1134, 3486, 120, 126, 1332,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3882, 0, 3, 1152, 3504, 138, 144, 1350,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3918, 0, 3, 1161, 3522, 144, 150, 1368,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3954, 0, 3, 1170, 3540, 150, 156, 1386,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3990, 0, 3, 1179, 3558, 156, 162, 1404,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4026, 0, 3, 1188, 3576, 162, 168, 1422,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4062, 0, 3, 1197, 3594, 168, 174, 1440,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4098, 0, 3, 1206, 3612, 174, 180, 1458,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4134, 0, 3, 1242, 3702, 192, 202, 1506,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4194, 0, 3, 1260, 3738, 202, 212, 1536,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4254, 0, 3, 1278, 3774, 212, 222, 1566,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4314, 0, 3, 1296, 3810, 222, 232, 1596,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4374, 0, 3, 1314, 3846, 232, 242, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4434, 0, 3, 1368, 3954, 262, 272, 1686,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4494, 0, 3, 1386, 3990, 272, 282, 1716,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4554, 0, 3, 1404, 4026, 282, 292, 1746,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4614, 0, 3, 1422, 4062, 292, 302, 1776,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4674, 0, 3, 1440, 4098, 302, 312, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4734, 0, 3, 3630, 3666, 1476, 4134, 332,
                                                 347, 1836, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4824, 0, 3, 3666, 3702, 1506, 4194, 347,
                                                 362, 1881, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4914, 0, 3, 3702, 3738, 1536, 4254, 362,
                                                 377, 1926, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5004, 0, 3, 3738, 3774, 1566, 4314, 377,
                                                 392, 1971, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5094, 0, 3, 3774, 3810, 1596, 4374, 392,
                                                 407, 2016, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5184, 0, 3, 3882, 3918, 1656, 4434, 437,
                                                 452, 2061, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5274, 0, 3, 3918, 3954, 1686, 4494, 452,
                                                 467, 2106, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5364, 0, 3, 3954, 3990, 1716, 4554, 467,
                                                 482, 2151, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5454, 0, 3, 3990, 4026, 1746, 4614, 482,
                                                 497, 2196, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5544, 0, 3, 4026, 4062, 1776, 4674, 497,
                                                 512, 2241, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5634, 0, 3, 4134, 4194, 1881, 4914, 542,
                                                 563, 2349, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5760, 0, 3, 4194, 4254, 1926, 5004, 563,
                                                 584, 2412, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5886, 0, 3, 4254, 4314, 1971, 5094, 584,
                                                 605, 2475, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6012, 0, 3, 4434, 4494, 2106, 5364, 647,
                                                 668, 2601, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6138, 0, 3, 4494, 4554, 2151, 5454, 668,
                                                 689, 2664, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6264, 0, 3, 4554, 4614, 2196, 5544, 689,
                                                 710, 2727, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6390, 0, 3, 4734, 4824, 2286, 5634, 752,
                                                 780, 2790, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6558, 0, 3, 4824, 4914, 2349, 5760, 780,
                                                 808, 2874, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6726, 0, 3, 4914, 5004, 2412, 5886, 808,
                                                 836, 2958, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6894, 0, 3, 5184, 5274, 2538, 6012, 892,
                                                 920, 3042, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7062, 0, 3, 5274, 5364, 2601, 6138, 920,
                                                 948, 3126, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7230, 0, 3, 5364, 5454, 2664, 6264, 948,
                                                 976, 3210, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7398, 3, 1032, 1035, 3300, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7408, 3, 1035, 1038, 3306, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7418, 3, 1038, 1041, 3312, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7428, 3, 1041, 1044, 3318, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7438, 3, 1044, 1047, 3324, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7448, 3, 1047, 1050, 3330, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7458, 3, 1056, 1059, 3342, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7468, 3, 1059, 1062, 3348, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7478, 3, 1062, 1065, 3354, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7488, 3, 1065, 1068, 3360, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7498, 3, 1068, 1071, 3366, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7508, 3, 1071, 1074, 3372, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 7518, 0, 3, 3294, 7398, 3396, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7548, 0, 3, 3300, 7408, 3414, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7578, 0, 3, 3306, 7418, 3432, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7608, 0, 3, 3312, 7428, 3450, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7638, 0, 3, 3318, 7438, 3468, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7668, 0, 3, 3324, 7448, 3486, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7698, 0, 3, 3336, 7458, 3522, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7728, 0, 3, 3342, 7468, 3540, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7758, 0, 3, 3348, 7478, 3558, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7788, 0, 3, 3354, 7488, 3576, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7818, 0, 3, 3360, 7498, 3594, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7848, 0, 3, 3366, 7508, 3612, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 7878, 0, 3, 3396, 7548, 1224, 1242,
                                                 3702, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7938, 0, 3, 3414, 7578, 1242, 1260,
                                                 3738, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7998, 0, 3, 3432, 7608, 1260, 1278,
                                                 3774, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8058, 0, 3, 3450, 7638, 1278, 1296,
                                                 3810, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8118, 0, 3, 3468, 7668, 1296, 1314,
                                                 3846, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8178, 0, 3, 3522, 7728, 1350, 1368,
                                                 3954, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8238, 0, 3, 3540, 7758, 1368, 1386,
                                                 3990, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8298, 0, 3, 3558, 7788, 1386, 1404,
                                                 4026, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8358, 0, 3, 3576, 7818, 1404, 1422,
                                                 4062, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8418, 0, 3, 3594, 7848, 1422, 1440,
                                                 4098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8478, 0, 3, 3702, 7938, 1476, 1506,
                                                 4194, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8578, 0, 3, 3738, 7998, 1506, 1536,
                                                 4254, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8678, 0, 3, 3774, 8058, 1536, 1566,
                                                 4314, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8778, 0, 3, 3810, 8118, 1566, 1596,
                                                 4374, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8878, 0, 3, 3954, 8238, 1656, 1686,
                                                 4494, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8978, 0, 3, 3990, 8298, 1686, 1716,
                                                 4554, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9078, 0, 3, 4026, 8358, 1716, 1746,
                                                 4614, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9178, 0, 3, 4062, 8418, 1746, 1776,
                                                 4674, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9278, 0, 3, 7878, 7938, 4194, 8578,
                                                 1836, 1881, 4914, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9428, 0, 3, 7938, 7998, 4254, 8678,
                                                 1881, 1926, 5004, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9578, 0, 3, 7998, 8058, 4314, 8778,
                                                 1926, 1971, 5094, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9728, 0, 3, 8178, 8238, 4494, 8978,
                                                 2061, 2106, 5364, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9878, 0, 3, 8238, 8298, 4554, 9078,
                                                 2106, 2151, 5454, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10028, 0, 3, 8298, 8358, 4614, 9178,
                                                 2151, 2196, 5544, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10178, 0, 3, 8478, 8578, 4914, 9428,
                                                 2286, 2349, 5760, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10388, 0, 3, 8578, 8678, 5004, 9578,
                                                 2349, 2412, 5886, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10598, 0, 3, 8878, 8978, 5364, 9878,
                                                 2538, 2601, 6138, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10808, 0, 3, 8978, 9078, 5454, 10028,
                                                 2601, 2664, 6264, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11018, 0, 3, 9278, 9428, 5760, 10388,
                                                 2790, 2874, 6726, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11298, 0, 3, 9728, 9878, 6138, 10808,
                                                 3042, 3126, 7230, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11578, 3, 3294, 3300, 7408, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11593, 3, 3300, 3306, 7418, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11608, 3, 3306, 3312, 7428, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11623, 3, 3312, 3318, 7438, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11638, 3, 3318, 3324, 7448, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11653, 3, 3336, 3342, 7468, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11668, 3, 3342, 3348, 7478, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11683, 3, 3348, 3354, 7488, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11698, 3, 3354, 3360, 7498, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11713, 3, 3360, 3366, 7508, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 11728, 0, 3, 7398, 11578, 7548, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11773, 0, 3, 7408, 11593, 7578, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11818, 0, 3, 7418, 11608, 7608, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11863, 0, 3, 7428, 11623, 7638, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11908, 0, 3, 7438, 11638, 7668, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11953, 0, 3, 7458, 11653, 7728, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11998, 0, 3, 7468, 11668, 7758, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12043, 0, 3, 7478, 11683, 7788, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12088, 0, 3, 7488, 11698, 7818, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12133, 0, 3, 7498, 11713, 7848, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 12178, 0, 3, 7518, 11728, 3630, 3666,
                                                 7878, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12268, 0, 3, 7548, 11773, 3666, 3702,
                                                 7938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12358, 0, 3, 7578, 11818, 3702, 3738,
                                                 7998, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12448, 0, 3, 7608, 11863, 3738, 3774,
                                                 8058, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12538, 0, 3, 7638, 11908, 3774, 3810,
                                                 8118, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12628, 0, 3, 7698, 11953, 3882, 3918,
                                                 8178, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12718, 0, 3, 7728, 11998, 3918, 3954,
                                                 8238, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12808, 0, 3, 7758, 12043, 3954, 3990,
                                                 8298, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12898, 0, 3, 7788, 12088, 3990, 4026,
                                                 8358, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12988, 0, 3, 7818, 12133, 4026, 4062,
                                                 8418, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13078, 0, 3, 7938, 12358, 4134, 4194,
                                                 8578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13228, 0, 3, 7998, 12448, 4194, 4254,
                                                 8678, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13378, 0, 3, 8058, 12538, 4254, 4314,
                                                 8778, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13528, 0, 3, 8238, 12808, 4434, 4494,
                                                 8978, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13678, 0, 3, 8298, 12898, 4494, 4554,
                                                 9078, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13828, 0, 3, 8358, 12988, 4554, 4614,
                                                 9178, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13978, 0, 3, 12178, 12268, 8478, 13078,
                                                 4734, 4824, 9278, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14203, 0, 3, 12268, 12358, 8578, 13228,
                                                 4824, 4914, 9428, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14428, 0, 3, 12358, 12448, 8678, 13378,
                                                 4914, 5004, 9578, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14653, 0, 3, 12628, 12718, 8878, 13528,
                                                 5184, 5274, 9728, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14878, 0, 3, 12718, 12808, 8978, 13678,
                                                 5274, 5364, 9878, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15103, 0, 3, 12808, 12898, 9078, 13828,
                                                 5364, 5454, 10028, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 15328, 0, 3, 13078, 13228, 9428, 14428,
                                                 5634, 5760, 10388, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 15643, 0, 3, 13528, 13678, 9878, 15103,
                                                 6012, 6138, 10808, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 15958, 0, 3, 13978, 14203, 10178, 15328,
                                                 6390, 6558, 11018, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 16378, 0, 3, 14653, 14878, 10598, 15643,
                                                 6894, 7062, 11298, ncols, alpha, beta, p);

            simdgeo::geom_h_x(buffer, 16798, 14653, 16378, 1, 15, ncols, alpha);

            simdgeo::geom_h_y(buffer, 17113, 14653, 16378, 1, 15, ncols, alpha);

            simdgeo::geom_h_z(buffer, 17428, 14653, 16378, 1, 15, ncols, alpha);

            simdgeo::geom_h_x(buffer, 17743, 13978, 15958, 1, 15, ncols, alpha);

            simdgeo::geom_h_y(buffer, 18058, 13978, 15958, 1, 15, ncols, alpha);

            simdgeo::geom_h_z(buffer, 18373, 13978, 15958, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 18688, 17743, 945, ncols);

            simdfunc::contract_primitives(buffer, 19633, 16798, 945, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 20578, 19633, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 20578, 9, nmax);

    simdtrf::transform_g_inner(buffer, 20578, 19948, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 99 * nvalues, nvalues, buffer, 20578, 9, nmax);

    simdtrf::transform_g_inner(buffer, 20578, 20263, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 198 * nvalues, nvalues, buffer, 20578, 9, nmax);

    simdtrf::transform_g_inner(buffer, 20578, 18688, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 297 * nvalues, nvalues, buffer, 20578, 9, nmax);

    simdtrf::transform_g_inner(buffer, 20578, 19003, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 396 * nvalues, nvalues, buffer, 20578, 9, nmax);

    simdtrf::transform_g_inner(buffer, 20578, 19318, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 495 * nvalues, nvalues, buffer, 20578, 9, nmax);
}

}  // namespace simdt2ceri
