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


#include "SimdElectronRepulsionRecKK.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
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
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 56338, nvalues);

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
                                            10, 11, 12, 13, 14}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 20, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 16, 17, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 17, 18, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 18, 19, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 21, 24, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 39, 42, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 42, 45, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 45, 48, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 48, 51, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 51, 54, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 242, 0, 60, 66, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 257, 0, 66, 72, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 272, 0, 72, 78, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 287, 0, 78, 84, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 302, 0, 84, 90, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 317, 0, 90, 96, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 332, 0, 96, 102, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 347, 0, 102, 108, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 362, 0, 108, 114, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 377, 0, 114, 120, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 392, 0, 132, 142, 257, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 413, 0, 142, 152, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 434, 0, 152, 162, 287, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 455, 0, 162, 172, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 476, 0, 172, 182, 317, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 497, 0, 182, 192, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 518, 0, 192, 202, 347, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 539, 0, 202, 212, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 560, 0, 212, 222, 377, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 242, 257, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 257, 272, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 272, 287, 455, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 287, 302, 476, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 302, 317, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 721, 0, 317, 332, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 749, 0, 332, 347, 539, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 777, 0, 347, 362, 560, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 805, 0, 392, 413, 609, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 841, 0, 413, 434, 637, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 877, 0, 434, 455, 665, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 913, 0, 455, 476, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 949, 0, 476, 497, 721, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 985, 0, 497, 518, 749, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1021, 0, 518, 539, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1057, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1060, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1063, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1066, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1069, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1072, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1075, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1078, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1081, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1084, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1087, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1090, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1099, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1108, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1117, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1126, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1135, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1144, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1153, 3, 16, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1162, 3, 17, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1171, 3, 18, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1180, 3, 19, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1189, 0, 3, 24, 1090, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1207, 0, 3, 27, 1099, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1225, 0, 3, 30, 1108, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1243, 0, 3, 33, 1117, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1261, 0, 3, 36, 1126, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1279, 0, 3, 39, 1135, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1297, 0, 3, 42, 1144, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1315, 0, 3, 45, 1153, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1333, 0, 3, 48, 1162, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1351, 0, 3, 51, 1171, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1369, 0, 3, 54, 1180, 126, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1387, 0, 3, 60, 1189, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1417, 0, 3, 66, 1207, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1447, 0, 3, 72, 1225, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1477, 0, 3, 78, 1243, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1507, 0, 3, 84, 1261, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1537, 0, 3, 90, 1279, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1567, 0, 3, 96, 1297, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1597, 0, 3, 102, 1315, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1627, 0, 3, 108, 1333, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1657, 0, 3, 114, 1351, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1687, 0, 3, 120, 1369, 232, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1717, 0, 3, 142, 1447, 257, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1762, 0, 3, 152, 1477, 272, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1807, 0, 3, 162, 1507, 287, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1852, 0, 3, 172, 1537, 302, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1897, 0, 3, 182, 1567, 317, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1942, 0, 3, 192, 1597, 332, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1987, 0, 3, 202, 1627, 347, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2032, 0, 3, 212, 1657, 362, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2077, 0, 3, 222, 1687, 377, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2122, 0, 3, 242, 1717, 392, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2185, 0, 3, 257, 1762, 413, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2248, 0, 3, 272, 1807, 434, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2311, 0, 3, 287, 1852, 455, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2374, 0, 3, 302, 1897, 476, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2437, 0, 3, 317, 1942, 497, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2500, 0, 3, 332, 1987, 518, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2563, 0, 3, 347, 2032, 539, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2626, 0, 3, 362, 2077, 560, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2689, 0, 3, 413, 2248, 609, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2773, 0, 3, 434, 2311, 637, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2857, 0, 3, 455, 2374, 665, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2941, 0, 3, 476, 2437, 693, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3025, 0, 3, 497, 2500, 721, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3109, 0, 3, 518, 2563, 749, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3193, 0, 3, 539, 2626, 777, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3277, 0, 3, 581, 2689, 805, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3385, 0, 3, 609, 2773, 841, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3493, 0, 3, 637, 2857, 877, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3601, 0, 3, 665, 2941, 913, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3709, 0, 3, 693, 3025, 949, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3817, 0, 3, 721, 3109, 985, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3925, 0, 3, 749, 3193, 1021, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4033, 3, 9, 10, 1060, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4039, 3, 10, 11, 1063, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4045, 3, 11, 12, 1066, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4051, 3, 12, 13, 1069, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4057, 3, 13, 14, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4063, 3, 14, 15, 1075, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4069, 3, 15, 16, 1078, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4075, 3, 16, 17, 1081, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4081, 3, 17, 18, 1084, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4087, 3, 18, 19, 1087, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4093, 0, 3, 1057, 4033, 1099, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4111, 0, 3, 1060, 4039, 1108, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4129, 0, 3, 1063, 4045, 1117, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4147, 0, 3, 1066, 4051, 1126, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4165, 0, 3, 1069, 4057, 1135, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4183, 0, 3, 1072, 4063, 1144, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4201, 0, 3, 1075, 4069, 1153, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4219, 0, 3, 1078, 4075, 1162, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4237, 0, 3, 1081, 4081, 1171, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4255, 0, 3, 1084, 4087, 1180, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4273, 0, 3, 1090, 4093, 60, 66, 1207,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4309, 0, 3, 1099, 4111, 66, 72, 1225,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4345, 0, 3, 1108, 4129, 72, 78, 1243,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4381, 0, 3, 1117, 4147, 78, 84, 1261,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4417, 0, 3, 1126, 4165, 84, 90, 1279,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4453, 0, 3, 1135, 4183, 90, 96, 1297,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4489, 0, 3, 1144, 4201, 96, 102, 1315,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4525, 0, 3, 1153, 4219, 102, 108, 1333,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4561, 0, 3, 1162, 4237, 108, 114, 1351,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4597, 0, 3, 1171, 4255, 114, 120, 1369,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4633, 0, 3, 1207, 4309, 132, 142, 1447,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4693, 0, 3, 1225, 4345, 142, 152, 1477,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4753, 0, 3, 1243, 4381, 152, 162, 1507,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4813, 0, 3, 1261, 4417, 162, 172, 1537,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4873, 0, 3, 1279, 4453, 172, 182, 1567,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4933, 0, 3, 1297, 4489, 182, 192, 1597,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4993, 0, 3, 1315, 4525, 192, 202, 1627,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5053, 0, 3, 1333, 4561, 202, 212, 1657,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5113, 0, 3, 1351, 4597, 212, 222, 1687,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5173, 0, 3, 4273, 4309, 1447, 4693, 242,
                                                 257, 1762, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5263, 0, 3, 4309, 4345, 1477, 4753, 257,
                                                 272, 1807, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5353, 0, 3, 4345, 4381, 1507, 4813, 272,
                                                 287, 1852, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5443, 0, 3, 4381, 4417, 1537, 4873, 287,
                                                 302, 1897, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5533, 0, 3, 4417, 4453, 1567, 4933, 302,
                                                 317, 1942, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5623, 0, 3, 4453, 4489, 1597, 4993, 317,
                                                 332, 1987, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5713, 0, 3, 4489, 4525, 1627, 5053, 332,
                                                 347, 2032, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5803, 0, 3, 4525, 4561, 1657, 5113, 347,
                                                 362, 2077, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5893, 0, 3, 4633, 4693, 1762, 5263, 392,
                                                 413, 2248, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6019, 0, 3, 4693, 4753, 1807, 5353, 413,
                                                 434, 2311, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6145, 0, 3, 4753, 4813, 1852, 5443, 434,
                                                 455, 2374, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6271, 0, 3, 4813, 4873, 1897, 5533, 455,
                                                 476, 2437, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6397, 0, 3, 4873, 4933, 1942, 5623, 476,
                                                 497, 2500, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6523, 0, 3, 4933, 4993, 1987, 5713, 497,
                                                 518, 2563, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6649, 0, 3, 4993, 5053, 2032, 5803, 518,
                                                 539, 2626, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6775, 0, 3, 5173, 5263, 2248, 6019, 581,
                                                 609, 2773, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6943, 0, 3, 5263, 5353, 2311, 6145, 609,
                                                 637, 2857, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7111, 0, 3, 5353, 5443, 2374, 6271, 637,
                                                 665, 2941, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7279, 0, 3, 5443, 5533, 2437, 6397, 665,
                                                 693, 3025, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7447, 0, 3, 5533, 5623, 2500, 6523, 693,
                                                 721, 3109, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7615, 0, 3, 5623, 5713, 2563, 6649, 721,
                                                 749, 3193, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7783, 0, 3, 5893, 6019, 2773, 6943, 805,
                                                 841, 3493, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7999, 0, 3, 6019, 6145, 2857, 7111, 841,
                                                 877, 3601, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8215, 0, 3, 6145, 6271, 2941, 7279, 877,
                                                 913, 3709, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8431, 0, 3, 6271, 6397, 3025, 7447, 913,
                                                 949, 3817, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8647, 0, 3, 6397, 6523, 3109, 7615, 949,
                                                 985, 3925, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8863, 3, 1057, 1060, 4039, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8873, 3, 1060, 1063, 4045, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8883, 3, 1063, 1066, 4051, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8893, 3, 1066, 1069, 4057, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8903, 3, 1069, 1072, 4063, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8913, 3, 1072, 1075, 4069, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8923, 3, 1075, 1078, 4075, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8933, 3, 1078, 1081, 4081, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8943, 3, 1081, 1084, 4087, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8953, 0, 3, 4033, 8863, 4111, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8983, 0, 3, 4039, 8873, 4129, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9013, 0, 3, 4045, 8883, 4147, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9043, 0, 3, 4051, 8893, 4165, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9073, 0, 3, 4057, 8903, 4183, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9103, 0, 3, 4063, 8913, 4201, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9133, 0, 3, 4069, 8923, 4219, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9163, 0, 3, 4075, 8933, 4237, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9193, 0, 3, 4081, 8943, 4255, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 9223, 0, 3, 4093, 8953, 1189, 1207,
                                                 4309, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9283, 0, 3, 4111, 8983, 1207, 1225,
                                                 4345, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9343, 0, 3, 4129, 9013, 1225, 1243,
                                                 4381, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9403, 0, 3, 4147, 9043, 1243, 1261,
                                                 4417, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9463, 0, 3, 4165, 9073, 1261, 1279,
                                                 4453, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9523, 0, 3, 4183, 9103, 1279, 1297,
                                                 4489, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9583, 0, 3, 4201, 9133, 1297, 1315,
                                                 4525, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9643, 0, 3, 4219, 9163, 1315, 1333,
                                                 4561, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9703, 0, 3, 4237, 9193, 1333, 1351,
                                                 4597, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9763, 0, 3, 4273, 9223, 1387, 1417,
                                                 4633, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9863, 0, 3, 4309, 9283, 1417, 1447,
                                                 4693, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9963, 0, 3, 4345, 9343, 1447, 1477,
                                                 4753, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10063, 0, 3, 4381, 9403, 1477, 1507,
                                                 4813, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10163, 0, 3, 4417, 9463, 1507, 1537,
                                                 4873, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10263, 0, 3, 4453, 9523, 1537, 1567,
                                                 4933, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10363, 0, 3, 4489, 9583, 1567, 1597,
                                                 4993, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10463, 0, 3, 4525, 9643, 1597, 1627,
                                                 5053, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10563, 0, 3, 4561, 9703, 1627, 1657,
                                                 5113, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10663, 0, 3, 9223, 9283, 4693, 9963,
                                                 1717, 1762, 5263, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10813, 0, 3, 9283, 9343, 4753, 10063,
                                                 1762, 1807, 5353, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10963, 0, 3, 9343, 9403, 4813, 10163,
                                                 1807, 1852, 5443, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11113, 0, 3, 9403, 9463, 4873, 10263,
                                                 1852, 1897, 5533, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11263, 0, 3, 9463, 9523, 4933, 10363,
                                                 1897, 1942, 5623, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11413, 0, 3, 9523, 9583, 4993, 10463,
                                                 1942, 1987, 5713, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11563, 0, 3, 9583, 9643, 5053, 10563,
                                                 1987, 2032, 5803, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11713, 0, 3, 9763, 9863, 5173, 10663,
                                                 2122, 2185, 5893, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11923, 0, 3, 9863, 9963, 5263, 10813,
                                                 2185, 2248, 6019, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12133, 0, 3, 9963, 10063, 5353, 10963,
                                                 2248, 2311, 6145, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12343, 0, 3, 10063, 10163, 5443, 11113,
                                                 2311, 2374, 6271, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12553, 0, 3, 10163, 10263, 5533, 11263,
                                                 2374, 2437, 6397, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12763, 0, 3, 10263, 10363, 5623, 11413,
                                                 2437, 2500, 6523, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12973, 0, 3, 10363, 10463, 5713, 11563,
                                                 2500, 2563, 6649, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13183, 0, 3, 10663, 10813, 6019, 12133,
                                                 2689, 2773, 6943, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13463, 0, 3, 10813, 10963, 6145, 12343,
                                                 2773, 2857, 7111, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 13743, 0, 3, 10963, 11113, 6271, 12553,
                                                 2857, 2941, 7279, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14023, 0, 3, 11113, 11263, 6397, 12763,
                                                 2941, 3025, 7447, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14303, 0, 3, 11263, 11413, 6523, 12973,
                                                 3025, 3109, 7615, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 14583, 0, 3, 11713, 11923, 6775, 13183,
                                                 3277, 3385, 7783, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 14943, 0, 3, 11923, 12133, 6943, 13463,
                                                 3385, 3493, 7999, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 15303, 0, 3, 12133, 12343, 7111, 13743,
                                                 3493, 3601, 8215, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 15663, 0, 3, 12343, 12553, 7279, 14023,
                                                 3601, 3709, 8431, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16023, 0, 3, 12553, 12763, 7447, 14303,
                                                 3709, 3817, 8647, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16383, 3, 4033, 4039, 8873, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16398, 3, 4039, 4045, 8883, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16413, 3, 4045, 4051, 8893, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16428, 3, 4051, 4057, 8903, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16443, 3, 4057, 4063, 8913, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16458, 3, 4063, 4069, 8923, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16473, 3, 4069, 4075, 8933, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16488, 3, 4075, 4081, 8943, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 16503, 0, 3, 8863, 16383, 8983, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16548, 0, 3, 8873, 16398, 9013, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16593, 0, 3, 8883, 16413, 9043, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16638, 0, 3, 8893, 16428, 9073, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16683, 0, 3, 8903, 16443, 9103, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16728, 0, 3, 8913, 16458, 9133, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16773, 0, 3, 8923, 16473, 9163, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16818, 0, 3, 8933, 16488, 9193, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 16863, 0, 3, 8953, 16503, 4273, 4309,
                                                 9283, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16953, 0, 3, 8983, 16548, 4309, 4345,
                                                 9343, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17043, 0, 3, 9013, 16593, 4345, 4381,
                                                 9403, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17133, 0, 3, 9043, 16638, 4381, 4417,
                                                 9463, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17223, 0, 3, 9073, 16683, 4417, 4453,
                                                 9523, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17313, 0, 3, 9103, 16728, 4453, 4489,
                                                 9583, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17403, 0, 3, 9133, 16773, 4489, 4525,
                                                 9643, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17493, 0, 3, 9163, 16818, 4525, 4561,
                                                 9703, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17583, 0, 3, 9283, 16953, 4633, 4693,
                                                 9963, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17733, 0, 3, 9343, 17043, 4693, 4753,
                                                 10063, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17883, 0, 3, 9403, 17133, 4753, 4813,
                                                 10163, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18033, 0, 3, 9463, 17223, 4813, 4873,
                                                 10263, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18183, 0, 3, 9523, 17313, 4873, 4933,
                                                 10363, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18333, 0, 3, 9583, 17403, 4933, 4993,
                                                 10463, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18483, 0, 3, 9643, 17493, 4993, 5053,
                                                 10563, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18633, 0, 3, 16863, 16953, 9963, 17733,
                                                 5173, 5263, 10813, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18858, 0, 3, 16953, 17043, 10063, 17883,
                                                 5263, 5353, 10963, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19083, 0, 3, 17043, 17133, 10163, 18033,
                                                 5353, 5443, 11113, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19308, 0, 3, 17133, 17223, 10263, 18183,
                                                 5443, 5533, 11263, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19533, 0, 3, 17223, 17313, 10363, 18333,
                                                 5533, 5623, 11413, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 19758, 0, 3, 17313, 17403, 10463, 18483,
                                                 5623, 5713, 11563, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 19983, 0, 3, 17583, 17733, 10813, 18858,
                                                 5893, 6019, 12133, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20298, 0, 3, 17733, 17883, 10963, 19083,
                                                 6019, 6145, 12343, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20613, 0, 3, 17883, 18033, 11113, 19308,
                                                 6145, 6271, 12553, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20928, 0, 3, 18033, 18183, 11263, 19533,
                                                 6271, 6397, 12763, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 21243, 0, 3, 18183, 18333, 11413, 19758,
                                                 6397, 6523, 12973, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 21558, 0, 3, 18633, 18858, 12133, 20298,
                                                 6775, 6943, 13463, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 21978, 0, 3, 18858, 19083, 12343, 20613,
                                                 6943, 7111, 13743, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 22398, 0, 3, 19083, 19308, 12553, 20928,
                                                 7111, 7279, 14023, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 22818, 0, 3, 19308, 19533, 12763, 21243,
                                                 7279, 7447, 14303, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 23238, 0, 3, 19983, 20298, 13463, 21978,
                                                 7783, 7999, 15303, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 23778, 0, 3, 20298, 20613, 13743, 22398,
                                                 7999, 8215, 15663, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 24318, 0, 3, 20613, 20928, 14023, 22818,
                                                 8215, 8431, 16023, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24858, 3, 8863, 8873, 16398, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24879, 3, 8873, 8883, 16413, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24900, 3, 8883, 8893, 16428, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24921, 3, 8893, 8903, 16443, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24942, 3, 8903, 8913, 16458, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24963, 3, 8913, 8923, 16473, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 24984, 3, 8923, 8933, 16488, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 25005, 0, 3, 16383, 24858, 16548, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25068, 0, 3, 16398, 24879, 16593, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25131, 0, 3, 16413, 24900, 16638, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25194, 0, 3, 16428, 24921, 16683, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25257, 0, 3, 16443, 24942, 16728, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25320, 0, 3, 16458, 24963, 16773, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25383, 0, 3, 16473, 24984, 16818, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 25446, 0, 3, 16503, 25005, 9223, 9283,
                                                 16953, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25572, 0, 3, 16548, 25068, 9283, 9343,
                                                 17043, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25698, 0, 3, 16593, 25131, 9343, 9403,
                                                 17133, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25824, 0, 3, 16638, 25194, 9403, 9463,
                                                 17223, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 25950, 0, 3, 16683, 25257, 9463, 9523,
                                                 17313, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 26076, 0, 3, 16728, 25320, 9523, 9583,
                                                 17403, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 26202, 0, 3, 16773, 25383, 9583, 9643,
                                                 17493, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26328, 0, 3, 16863, 25446, 9763, 9863,
                                                 17583, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26538, 0, 3, 16953, 25572, 9863, 9963,
                                                 17733, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26748, 0, 3, 17043, 25698, 9963, 10063,
                                                 17883, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 26958, 0, 3, 17133, 25824, 10063, 10163,
                                                 18033, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 27168, 0, 3, 17223, 25950, 10163, 10263,
                                                 18183, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 27378, 0, 3, 17313, 26076, 10263, 10363,
                                                 18333, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 27588, 0, 3, 17403, 26202, 10363, 10463,
                                                 18483, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 27798, 0, 3, 25446, 25572, 17733, 26748,
                                                 10663, 10813, 18858, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 28113, 0, 3, 25572, 25698, 17883, 26958,
                                                 10813, 10963, 19083, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 28428, 0, 3, 25698, 25824, 18033, 27168,
                                                 10963, 11113, 19308, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 28743, 0, 3, 25824, 25950, 18183, 27378,
                                                 11113, 11263, 19533, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 29058, 0, 3, 25950, 26076, 18333, 27588,
                                                 11263, 11413, 19758, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 29373, 0, 3, 26328, 26538, 18633, 27798,
                                                 11713, 11923, 19983, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 29814, 0, 3, 26538, 26748, 18858, 28113,
                                                 11923, 12133, 20298, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 30255, 0, 3, 26748, 26958, 19083, 28428,
                                                 12133, 12343, 20613, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 30696, 0, 3, 26958, 27168, 19308, 28743,
                                                 12343, 12553, 20928, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 31137, 0, 3, 27168, 27378, 19533, 29058,
                                                 12553, 12763, 21243, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 31578, 0, 3, 27798, 28113, 20298, 30255,
                                                 13183, 13463, 21978, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 32166, 0, 3, 28113, 28428, 20613, 30696,
                                                 13463, 13743, 22398, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 32754, 0, 3, 28428, 28743, 20928, 31137,
                                                 13743, 14023, 22818, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 33342, 0, 3, 29373, 29814, 21558, 31578,
                                                 14583, 14943, 23238, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 34098, 0, 3, 29814, 30255, 21978, 32166,
                                                 14943, 15303, 23778, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 34854, 0, 3, 30255, 30696, 22398, 32754,
                                                 15303, 15663, 24318, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35610, 3, 16383, 16398, 24879, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35638, 3, 16398, 16413, 24900, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35666, 3, 16413, 16428, 24921, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35694, 3, 16428, 16443, 24942, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35722, 3, 16443, 16458, 24963, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35750, 3, 16458, 16473, 24984, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 35778, 0, 3, 24858, 35610, 25068, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35862, 0, 3, 24879, 35638, 25131, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35946, 0, 3, 24900, 35666, 25194, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36030, 0, 3, 24921, 35694, 25257, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36114, 0, 3, 24942, 35722, 25320, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36198, 0, 3, 24963, 35750, 25383, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 36282, 0, 3, 25005, 35778, 16863, 16953,
                                                 25572, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36450, 0, 3, 25068, 35862, 16953, 17043,
                                                 25698, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36618, 0, 3, 25131, 35946, 17043, 17133,
                                                 25824, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36786, 0, 3, 25194, 36030, 17133, 17223,
                                                 25950, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36954, 0, 3, 25257, 36114, 17223, 17313,
                                                 26076, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37122, 0, 3, 25320, 36198, 17313, 17403,
                                                 26202, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 37290, 0, 3, 25572, 36450, 17583, 17733,
                                                 26748, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 37570, 0, 3, 25698, 36618, 17733, 17883,
                                                 26958, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 37850, 0, 3, 25824, 36786, 17883, 18033,
                                                 27168, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38130, 0, 3, 25950, 36954, 18033, 18183,
                                                 27378, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38410, 0, 3, 26076, 37122, 18183, 18333,
                                                 27588, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 38690, 0, 3, 36282, 36450, 26748, 37570,
                                                 18633, 18858, 28113, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 39110, 0, 3, 36450, 36618, 26958, 37850,
                                                 18858, 19083, 28428, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 39530, 0, 3, 36618, 36786, 27168, 38130,
                                                 19083, 19308, 28743, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 39950, 0, 3, 36786, 36954, 27378, 38410,
                                                 19308, 19533, 29058, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 40370, 0, 3, 37290, 37570, 28113, 39110,
                                                 19983, 20298, 30255, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 40958, 0, 3, 37570, 37850, 28428, 39530,
                                                 20298, 20613, 30696, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 41546, 0, 3, 37850, 38130, 28743, 39950,
                                                 20613, 20928, 31137, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 42134, 0, 3, 38690, 39110, 30255, 40958,
                                                 21558, 21978, 32166, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 42918, 0, 3, 39110, 39530, 30696, 41546,
                                                 21978, 22398, 32754, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 43702, 0, 3, 40370, 40958, 32166, 42918,
                                                 23238, 23778, 34854, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 44710, 3, 24858, 24879, 35638, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 44746, 3, 24879, 24900, 35666, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 44782, 3, 24900, 24921, 35694, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 44818, 3, 24921, 24942, 35722, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 44854, 3, 24942, 24963, 35750, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 44890, 0, 3, 35610, 44710, 35862, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 44998, 0, 3, 35638, 44746, 35946, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 45106, 0, 3, 35666, 44782, 36030, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 45214, 0, 3, 35694, 44818, 36114, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 45322, 0, 3, 35722, 44854, 36198, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 45430, 0, 3, 35778, 44890, 25446, 25572,
                                                 36450, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 45646, 0, 3, 35862, 44998, 25572, 25698,
                                                 36618, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 45862, 0, 3, 35946, 45106, 25698, 25824,
                                                 36786, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 46078, 0, 3, 36030, 45214, 25824, 25950,
                                                 36954, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 46294, 0, 3, 36114, 45322, 25950, 26076,
                                                 37122, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 46510, 0, 3, 36282, 45430, 26328, 26538,
                                                 37290, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 46870, 0, 3, 36450, 45646, 26538, 26748,
                                                 37570, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 47230, 0, 3, 36618, 45862, 26748, 26958,
                                                 37850, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 47590, 0, 3, 36786, 46078, 26958, 27168,
                                                 38130, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 47950, 0, 3, 36954, 46294, 27168, 27378,
                                                 38410, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 48310, 0, 3, 45430, 45646, 37570, 47230,
                                                 27798, 28113, 39110, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 48850, 0, 3, 45646, 45862, 37850, 47590,
                                                 28113, 28428, 39530, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 49390, 0, 3, 45862, 46078, 38130, 47950,
                                                 28428, 28743, 39950, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 49930, 0, 3, 46510, 46870, 38690, 48310,
                                                 29373, 29814, 40370, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 50686, 0, 3, 46870, 47230, 39110, 48850,
                                                 29814, 30255, 40958, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 51442, 0, 3, 47230, 47590, 39530, 49390,
                                                 30255, 30696, 41546, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 52198, 0, 3, 48310, 48850, 40958, 51442,
                                                 31578, 32166, 42918, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 53206, 0, 3, 49930, 50686, 42134, 52198,
                                                 33342, 34098, 43702, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 54502, 53206, 1296, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 55798, 54502, 36, nmax);

    simdtrf::transform_k_outer_tri(values, nvalues, buffer, 55798, nmax);
}

}  // namespace simdt2ceri
