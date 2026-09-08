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


#include "SimdElectronRepulsionRecKL.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIL.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKL.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
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
#include "SimdTransformK.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kl_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(82037, nvalues);

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
                                            10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 7, 8, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 8, 9, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 9, 10, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 85, 0, 10, 11, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 91, 0, 11, 12, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 97, 0, 12, 13, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 103, 0, 13, 14, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 109, 0, 14, 15, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 115, 0, 15, 16, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 121, 0, 16, 17, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 127, 0, 17, 18, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 133, 0, 18, 19, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 139, 0, 19, 20, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 22, 25, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 25, 28, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 28, 31, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 31, 34, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 34, 37, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 37, 40, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 40, 43, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 43, 46, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 46, 49, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 235, 0, 49, 52, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 245, 0, 52, 55, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 255, 0, 55, 58, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 265, 0, 58, 61, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 275, 0, 67, 73, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 290, 0, 73, 79, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 305, 0, 79, 85, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 320, 0, 85, 91, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 335, 0, 91, 97, 205, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 350, 0, 97, 103, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 365, 0, 103, 109, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 380, 0, 109, 115, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 395, 0, 115, 121, 245, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 410, 0, 121, 127, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 425, 0, 127, 133, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 440, 0, 145, 155, 275, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 461, 0, 155, 165, 290, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 482, 0, 165, 175, 305, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 503, 0, 175, 185, 320, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 524, 0, 185, 195, 335, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 545, 0, 195, 205, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 566, 0, 205, 215, 365, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 587, 0, 215, 225, 380, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 608, 0, 225, 235, 395, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 629, 0, 235, 245, 410, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 650, 0, 245, 255, 425, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 671, 0, 275, 290, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 699, 0, 290, 305, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 727, 0, 305, 320, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 755, 0, 320, 335, 545, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 783, 0, 335, 350, 566, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 811, 0, 350, 365, 587, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 839, 0, 365, 380, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 867, 0, 380, 395, 629, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 895, 0, 395, 410, 650, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 923, 0, 440, 461, 671, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 959, 0, 461, 482, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 995, 0, 482, 503, 727, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1031, 0, 503, 524, 755, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1067, 0, 524, 545, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1103, 0, 545, 566, 811, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1139, 0, 566, 587, 839, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1175, 0, 587, 608, 867, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1211, 0, 608, 629, 895, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1247, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1250, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1253, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1256, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1259, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1262, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1265, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1268, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1271, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1274, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1277, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1280, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1283, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1292, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1301, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1310, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1319, 3, 13, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1328, 3, 14, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1337, 3, 15, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1346, 3, 16, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1355, 3, 17, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1364, 3, 18, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1373, 3, 19, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1382, 3, 20, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1391, 0, 3, 28, 1283, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1409, 0, 3, 31, 1292, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1427, 0, 3, 34, 1301, 85, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1445, 0, 3, 37, 1310, 91, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1463, 0, 3, 40, 1319, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1481, 0, 3, 43, 1328, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1499, 0, 3, 46, 1337, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1517, 0, 3, 49, 1346, 115, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1535, 0, 3, 52, 1355, 121, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1553, 0, 3, 55, 1364, 127, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1571, 0, 3, 58, 1373, 133, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1589, 0, 3, 61, 1382, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1607, 0, 3, 73, 1409, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1637, 0, 3, 79, 1427, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1667, 0, 3, 85, 1445, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1697, 0, 3, 91, 1463, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1727, 0, 3, 97, 1481, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1757, 0, 3, 103, 1499, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1787, 0, 3, 109, 1517, 225, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1817, 0, 3, 115, 1535, 235, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1847, 0, 3, 121, 1553, 245, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1877, 0, 3, 127, 1571, 255, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1907, 0, 3, 133, 1589, 265, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1937, 0, 3, 165, 1637, 290, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1982, 0, 3, 175, 1667, 305, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2027, 0, 3, 185, 1697, 320, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2072, 0, 3, 195, 1727, 335, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2117, 0, 3, 205, 1757, 350, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2162, 0, 3, 215, 1787, 365, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2207, 0, 3, 225, 1817, 380, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2252, 0, 3, 235, 1847, 395, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2297, 0, 3, 245, 1877, 410, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2342, 0, 3, 255, 1907, 425, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2387, 0, 3, 290, 1982, 482, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2450, 0, 3, 305, 2027, 503, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2513, 0, 3, 320, 2072, 524, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2576, 0, 3, 335, 2117, 545, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2639, 0, 3, 350, 2162, 566, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2702, 0, 3, 365, 2207, 587, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2765, 0, 3, 380, 2252, 608, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2828, 0, 3, 395, 2297, 629, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2891, 0, 3, 410, 2342, 650, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2954, 0, 3, 482, 2450, 699, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3038, 0, 3, 503, 2513, 727, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3122, 0, 3, 524, 2576, 755, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3206, 0, 3, 545, 2639, 783, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3290, 0, 3, 566, 2702, 811, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3374, 0, 3, 587, 2765, 839, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3458, 0, 3, 608, 2828, 867, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3542, 0, 3, 629, 2891, 895, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3626, 0, 3, 699, 3038, 995, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3734, 0, 3, 727, 3122, 1031, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3842, 0, 3, 755, 3206, 1067, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3950, 0, 3, 783, 3290, 1103, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4058, 0, 3, 811, 3374, 1139, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4166, 0, 3, 839, 3458, 1175, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4274, 0, 3, 867, 3542, 1211, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4382, 3, 9, 10, 1250, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4388, 3, 10, 11, 1253, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4394, 3, 11, 12, 1256, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4400, 3, 12, 13, 1259, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4406, 3, 13, 14, 1262, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4412, 3, 14, 15, 1265, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4418, 3, 15, 16, 1268, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4424, 3, 16, 17, 1271, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4430, 3, 17, 18, 1274, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4436, 3, 18, 19, 1277, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4442, 3, 19, 20, 1280, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4448, 0, 3, 1247, 4382, 1292, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4466, 0, 3, 1250, 4388, 1301, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4484, 0, 3, 1253, 4394, 1310, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4502, 0, 3, 1256, 4400, 1319, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4520, 0, 3, 1259, 4406, 1328, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4538, 0, 3, 1262, 4412, 1337, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4556, 0, 3, 1265, 4418, 1346, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4574, 0, 3, 1268, 4424, 1355, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4592, 0, 3, 1271, 4430, 1364, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4610, 0, 3, 1274, 4436, 1373, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4628, 0, 3, 1277, 4442, 1382, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4646, 0, 3, 1283, 4448, 67, 73, 1409,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4682, 0, 3, 1292, 4466, 73, 79, 1427,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4718, 0, 3, 1301, 4484, 79, 85, 1445,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4754, 0, 3, 1310, 4502, 85, 91, 1463,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4790, 0, 3, 1319, 4520, 91, 97, 1481,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4826, 0, 3, 1328, 4538, 97, 103, 1499,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4862, 0, 3, 1337, 4556, 103, 109, 1517,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4898, 0, 3, 1346, 4574, 109, 115, 1535,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4934, 0, 3, 1355, 4592, 115, 121, 1553,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4970, 0, 3, 1364, 4610, 121, 127, 1571,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5006, 0, 3, 1373, 4628, 127, 133, 1589,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5042, 0, 3, 1391, 4646, 145, 155, 1607,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5102, 0, 3, 1409, 4682, 155, 165, 1637,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5162, 0, 3, 1427, 4718, 165, 175, 1667,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5222, 0, 3, 1445, 4754, 175, 185, 1697,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5282, 0, 3, 1463, 4790, 185, 195, 1727,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5342, 0, 3, 1481, 4826, 195, 205, 1757,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5402, 0, 3, 1499, 4862, 205, 215, 1787,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5462, 0, 3, 1517, 4898, 215, 225, 1817,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5522, 0, 3, 1535, 4934, 225, 235, 1847,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5582, 0, 3, 1553, 4970, 235, 245, 1877,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5642, 0, 3, 1571, 5006, 245, 255, 1907,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5702, 0, 3, 4646, 4682, 1637, 5162, 275,
                                                 290, 1982, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5792, 0, 3, 4682, 4718, 1667, 5222, 290,
                                                 305, 2027, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5882, 0, 3, 4718, 4754, 1697, 5282, 305,
                                                 320, 2072, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5972, 0, 3, 4754, 4790, 1727, 5342, 320,
                                                 335, 2117, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6062, 0, 3, 4790, 4826, 1757, 5402, 335,
                                                 350, 2162, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6152, 0, 3, 4826, 4862, 1787, 5462, 350,
                                                 365, 2207, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6242, 0, 3, 4862, 4898, 1817, 5522, 365,
                                                 380, 2252, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6332, 0, 3, 4898, 4934, 1847, 5582, 380,
                                                 395, 2297, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6422, 0, 3, 4934, 4970, 1877, 5642, 395,
                                                 410, 2342, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6512, 0, 3, 5042, 5102, 1937, 5702, 440,
                                                 461, 2387, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6638, 0, 3, 5102, 5162, 1982, 5792, 461,
                                                 482, 2450, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6764, 0, 3, 5162, 5222, 2027, 5882, 482,
                                                 503, 2513, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6890, 0, 3, 5222, 5282, 2072, 5972, 503,
                                                 524, 2576, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7016, 0, 3, 5282, 5342, 2117, 6062, 524,
                                                 545, 2639, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7142, 0, 3, 5342, 5402, 2162, 6152, 545,
                                                 566, 2702, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7268, 0, 3, 5402, 5462, 2207, 6242, 566,
                                                 587, 2765, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7394, 0, 3, 5462, 5522, 2252, 6332, 587,
                                                 608, 2828, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7520, 0, 3, 5522, 5582, 2297, 6422, 608,
                                                 629, 2891, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7646, 0, 3, 5702, 5792, 2450, 6764, 671,
                                                 699, 3038, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7814, 0, 3, 5792, 5882, 2513, 6890, 699,
                                                 727, 3122, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7982, 0, 3, 5882, 5972, 2576, 7016, 727,
                                                 755, 3206, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8150, 0, 3, 5972, 6062, 2639, 7142, 755,
                                                 783, 3290, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8318, 0, 3, 6062, 6152, 2702, 7268, 783,
                                                 811, 3374, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8486, 0, 3, 6152, 6242, 2765, 7394, 811,
                                                 839, 3458, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8654, 0, 3, 6242, 6332, 2828, 7520, 839,
                                                 867, 3542, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8822, 0, 3, 6512, 6638, 2954, 7646, 923,
                                                 959, 3626, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9038, 0, 3, 6638, 6764, 3038, 7814, 959,
                                                 995, 3734, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9254, 0, 3, 6764, 6890, 3122, 7982, 995,
                                                 1031, 3842, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9470, 0, 3, 6890, 7016, 3206, 8150,
                                                 1031, 1067, 3950, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9686, 0, 3, 7016, 7142, 3290, 8318,
                                                 1067, 1103, 4058, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9902, 0, 3, 7142, 7268, 3374, 8486,
                                                 1103, 1139, 4166, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10118, 0, 3, 7268, 7394, 3458, 8654,
                                                 1139, 1175, 4274, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10334, 3, 1247, 1250, 4388, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10344, 3, 1250, 1253, 4394, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10354, 3, 1253, 1256, 4400, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10364, 3, 1256, 1259, 4406, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10374, 3, 1259, 1262, 4412, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10384, 3, 1262, 1265, 4418, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10394, 3, 1265, 1268, 4424, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10404, 3, 1268, 1271, 4430, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10414, 3, 1271, 1274, 4436, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10424, 3, 1274, 1277, 4442, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 10434, 0, 3, 4382, 10334, 4466, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10464, 0, 3, 4388, 10344, 4484, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10494, 0, 3, 4394, 10354, 4502, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10524, 0, 3, 4400, 10364, 4520, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10554, 0, 3, 4406, 10374, 4538, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10584, 0, 3, 4412, 10384, 4556, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10614, 0, 3, 4418, 10394, 4574, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10644, 0, 3, 4424, 10404, 4592, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10674, 0, 3, 4430, 10414, 4610, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10704, 0, 3, 4436, 10424, 4628, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 10734, 0, 3, 4448, 10434, 1391, 1409,
                                                 4682, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10794, 0, 3, 4466, 10464, 1409, 1427,
                                                 4718, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10854, 0, 3, 4484, 10494, 1427, 1445,
                                                 4754, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10914, 0, 3, 4502, 10524, 1445, 1463,
                                                 4790, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10974, 0, 3, 4520, 10554, 1463, 1481,
                                                 4826, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11034, 0, 3, 4538, 10584, 1481, 1499,
                                                 4862, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11094, 0, 3, 4556, 10614, 1499, 1517,
                                                 4898, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11154, 0, 3, 4574, 10644, 1517, 1535,
                                                 4934, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11214, 0, 3, 4592, 10674, 1535, 1553,
                                                 4970, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11274, 0, 3, 4610, 10704, 1553, 1571,
                                                 5006, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11334, 0, 3, 4682, 10794, 1607, 1637,
                                                 5162, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11434, 0, 3, 4718, 10854, 1637, 1667,
                                                 5222, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11534, 0, 3, 4754, 10914, 1667, 1697,
                                                 5282, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11634, 0, 3, 4790, 10974, 1697, 1727,
                                                 5342, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11734, 0, 3, 4826, 11034, 1727, 1757,
                                                 5402, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11834, 0, 3, 4862, 11094, 1757, 1787,
                                                 5462, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11934, 0, 3, 4898, 11154, 1787, 1817,
                                                 5522, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12034, 0, 3, 4934, 11214, 1817, 1847,
                                                 5582, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12134, 0, 3, 4970, 11274, 1847, 1877,
                                                 5642, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12234, 0, 3, 10734, 10794, 5162, 11434,
                                                 1937, 1982, 5792, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12384, 0, 3, 10794, 10854, 5222, 11534,
                                                 1982, 2027, 5882, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12534, 0, 3, 10854, 10914, 5282, 11634,
                                                 2027, 2072, 5972, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12684, 0, 3, 10914, 10974, 5342, 11734,
                                                 2072, 2117, 6062, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12834, 0, 3, 10974, 11034, 5402, 11834,
                                                 2117, 2162, 6152, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12984, 0, 3, 11034, 11094, 5462, 11934,
                                                 2162, 2207, 6242, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13134, 0, 3, 11094, 11154, 5522, 12034,
                                                 2207, 2252, 6332, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13284, 0, 3, 11154, 11214, 5582, 12134,
                                                 2252, 2297, 6422, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13434, 0, 3, 11334, 11434, 5792, 12384,
                                                 2387, 2450, 6764, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13644, 0, 3, 11434, 11534, 5882, 12534,
                                                 2450, 2513, 6890, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13854, 0, 3, 11534, 11634, 5972, 12684,
                                                 2513, 2576, 7016, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14064, 0, 3, 11634, 11734, 6062, 12834,
                                                 2576, 2639, 7142, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14274, 0, 3, 11734, 11834, 6152, 12984,
                                                 2639, 2702, 7268, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14484, 0, 3, 11834, 11934, 6242, 13134,
                                                 2702, 2765, 7394, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14694, 0, 3, 11934, 12034, 6332, 13284,
                                                 2765, 2828, 7520, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14904, 0, 3, 12234, 12384, 6764, 13644,
                                                 2954, 3038, 7814, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15184, 0, 3, 12384, 12534, 6890, 13854,
                                                 3038, 3122, 7982, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15464, 0, 3, 12534, 12684, 7016, 14064,
                                                 3122, 3206, 8150, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15744, 0, 3, 12684, 12834, 7142, 14274,
                                                 3206, 3290, 8318, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16024, 0, 3, 12834, 12984, 7268, 14484,
                                                 3290, 3374, 8486, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16304, 0, 3, 12984, 13134, 7394, 14694,
                                                 3374, 3458, 8654, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16584, 0, 3, 13434, 13644, 7814, 15184,
                                                 3626, 3734, 9254, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16944, 0, 3, 13644, 13854, 7982, 15464,
                                                 3734, 3842, 9470, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17304, 0, 3, 13854, 14064, 8150, 15744,
                                                 3842, 3950, 9686, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17664, 0, 3, 14064, 14274, 8318, 16024,
                                                 3950, 4058, 9902, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 18024, 0, 3, 14274, 14484, 8486, 16304,
                                                 4058, 4166, 10118, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18384, 3, 4382, 4388, 10344, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18399, 3, 4388, 4394, 10354, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18414, 3, 4394, 4400, 10364, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18429, 3, 4400, 4406, 10374, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18444, 3, 4406, 4412, 10384, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18459, 3, 4412, 4418, 10394, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18474, 3, 4418, 4424, 10404, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18489, 3, 4424, 4430, 10414, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18504, 3, 4430, 4436, 10424, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 18519, 0, 3, 10334, 18384, 10464, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18564, 0, 3, 10344, 18399, 10494, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18609, 0, 3, 10354, 18414, 10524, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18654, 0, 3, 10364, 18429, 10554, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18699, 0, 3, 10374, 18444, 10584, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18744, 0, 3, 10384, 18459, 10614, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18789, 0, 3, 10394, 18474, 10644, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18834, 0, 3, 10404, 18489, 10674, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18879, 0, 3, 10414, 18504, 10704, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 18924, 0, 3, 10434, 18519, 4646, 4682,
                                                 10794, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19014, 0, 3, 10464, 18564, 4682, 4718,
                                                 10854, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19104, 0, 3, 10494, 18609, 4718, 4754,
                                                 10914, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19194, 0, 3, 10524, 18654, 4754, 4790,
                                                 10974, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19284, 0, 3, 10554, 18699, 4790, 4826,
                                                 11034, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19374, 0, 3, 10584, 18744, 4826, 4862,
                                                 11094, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19464, 0, 3, 10614, 18789, 4862, 4898,
                                                 11154, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19554, 0, 3, 10644, 18834, 4898, 4934,
                                                 11214, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19644, 0, 3, 10674, 18879, 4934, 4970,
                                                 11274, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19734, 0, 3, 10734, 18924, 5042, 5102,
                                                 11334, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19884, 0, 3, 10794, 19014, 5102, 5162,
                                                 11434, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20034, 0, 3, 10854, 19104, 5162, 5222,
                                                 11534, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20184, 0, 3, 10914, 19194, 5222, 5282,
                                                 11634, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20334, 0, 3, 10974, 19284, 5282, 5342,
                                                 11734, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20484, 0, 3, 11034, 19374, 5342, 5402,
                                                 11834, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20634, 0, 3, 11094, 19464, 5402, 5462,
                                                 11934, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20784, 0, 3, 11154, 19554, 5462, 5522,
                                                 12034, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20934, 0, 3, 11214, 19644, 5522, 5582,
                                                 12134, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21084, 0, 3, 18924, 19014, 11434, 20034,
                                                 5702, 5792, 12384, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21309, 0, 3, 19014, 19104, 11534, 20184,
                                                 5792, 5882, 12534, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21534, 0, 3, 19104, 19194, 11634, 20334,
                                                 5882, 5972, 12684, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21759, 0, 3, 19194, 19284, 11734, 20484,
                                                 5972, 6062, 12834, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21984, 0, 3, 19284, 19374, 11834, 20634,
                                                 6062, 6152, 12984, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22209, 0, 3, 19374, 19464, 11934, 20784,
                                                 6152, 6242, 13134, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22434, 0, 3, 19464, 19554, 12034, 20934,
                                                 6242, 6332, 13284, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22659, 0, 3, 19734, 19884, 12234, 21084,
                                                 6512, 6638, 13434, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22974, 0, 3, 19884, 20034, 12384, 21309,
                                                 6638, 6764, 13644, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23289, 0, 3, 20034, 20184, 12534, 21534,
                                                 6764, 6890, 13854, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23604, 0, 3, 20184, 20334, 12684, 21759,
                                                 6890, 7016, 14064, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23919, 0, 3, 20334, 20484, 12834, 21984,
                                                 7016, 7142, 14274, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24234, 0, 3, 20484, 20634, 12984, 22209,
                                                 7142, 7268, 14484, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24549, 0, 3, 20634, 20784, 13134, 22434,
                                                 7268, 7394, 14694, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 24864, 0, 3, 21084, 21309, 13644, 23289,
                                                 7646, 7814, 15184, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 25284, 0, 3, 21309, 21534, 13854, 23604,
                                                 7814, 7982, 15464, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 25704, 0, 3, 21534, 21759, 14064, 23919,
                                                 7982, 8150, 15744, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26124, 0, 3, 21759, 21984, 14274, 24234,
                                                 8150, 8318, 16024, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26544, 0, 3, 21984, 22209, 14484, 24549,
                                                 8318, 8486, 16304, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 26964, 0, 3, 22659, 22974, 14904, 24864,
                                                 8822, 9038, 16584, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 27504, 0, 3, 22974, 23289, 15184, 25284,
                                                 9038, 9254, 16944, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 28044, 0, 3, 23289, 23604, 15464, 25704,
                                                 9254, 9470, 17304, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 28584, 0, 3, 23604, 23919, 15744, 26124,
                                                 9470, 9686, 17664, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 29124, 0, 3, 23919, 24234, 16024, 26544,
                                                 9686, 9902, 18024, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29664, 3, 10334, 10344, 18399, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29685, 3, 10344, 10354, 18414, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29706, 3, 10354, 10364, 18429, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29727, 3, 10364, 10374, 18444, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29748, 3, 10374, 10384, 18459, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29769, 3, 10384, 10394, 18474, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29790, 3, 10394, 10404, 18489, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 29811, 3, 10404, 10414, 18504, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 29832, 0, 3, 18384, 29664, 18564, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29895, 0, 3, 18399, 29685, 18609, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29958, 0, 3, 18414, 29706, 18654, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 30021, 0, 3, 18429, 29727, 18699, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 30084, 0, 3, 18444, 29748, 18744, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 30147, 0, 3, 18459, 29769, 18789, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 30210, 0, 3, 18474, 29790, 18834, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 30273, 0, 3, 18489, 29811, 18879, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 30336, 0, 3, 18519, 29832, 10734, 10794,
                                                 19014, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30462, 0, 3, 18564, 29895, 10794, 10854,
                                                 19104, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30588, 0, 3, 18609, 29958, 10854, 10914,
                                                 19194, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30714, 0, 3, 18654, 30021, 10914, 10974,
                                                 19284, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30840, 0, 3, 18699, 30084, 10974, 11034,
                                                 19374, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30966, 0, 3, 18744, 30147, 11034, 11094,
                                                 19464, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 31092, 0, 3, 18789, 30210, 11094, 11154,
                                                 19554, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 31218, 0, 3, 18834, 30273, 11154, 11214,
                                                 19644, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31344, 0, 3, 19014, 30462, 11334, 11434,
                                                 20034, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31554, 0, 3, 19104, 30588, 11434, 11534,
                                                 20184, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31764, 0, 3, 19194, 30714, 11534, 11634,
                                                 20334, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31974, 0, 3, 19284, 30840, 11634, 11734,
                                                 20484, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32184, 0, 3, 19374, 30966, 11734, 11834,
                                                 20634, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32394, 0, 3, 19464, 31092, 11834, 11934,
                                                 20784, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32604, 0, 3, 19554, 31218, 11934, 12034,
                                                 20934, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 32814, 0, 3, 30336, 30462, 20034, 31554,
                                                 12234, 12384, 21309, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33129, 0, 3, 30462, 30588, 20184, 31764,
                                                 12384, 12534, 21534, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33444, 0, 3, 30588, 30714, 20334, 31974,
                                                 12534, 12684, 21759, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33759, 0, 3, 30714, 30840, 20484, 32184,
                                                 12684, 12834, 21984, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34074, 0, 3, 30840, 30966, 20634, 32394,
                                                 12834, 12984, 22209, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34389, 0, 3, 30966, 31092, 20784, 32604,
                                                 12984, 13134, 22434, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34704, 0, 3, 31344, 31554, 21309, 33129,
                                                 13434, 13644, 23289, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35145, 0, 3, 31554, 31764, 21534, 33444,
                                                 13644, 13854, 23604, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35586, 0, 3, 31764, 31974, 21759, 33759,
                                                 13854, 14064, 23919, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 36027, 0, 3, 31974, 32184, 21984, 34074,
                                                 14064, 14274, 24234, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 36468, 0, 3, 32184, 32394, 22209, 34389,
                                                 14274, 14484, 24549, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 36909, 0, 3, 32814, 33129, 23289, 35145,
                                                 14904, 15184, 25284, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 37497, 0, 3, 33129, 33444, 23604, 35586,
                                                 15184, 15464, 25704, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 38085, 0, 3, 33444, 33759, 23919, 36027,
                                                 15464, 15744, 26124, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 38673, 0, 3, 33759, 34074, 24234, 36468,
                                                 15744, 16024, 26544, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 39261, 0, 3, 34704, 35145, 25284, 37497,
                                                 16584, 16944, 28044, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 40017, 0, 3, 35145, 35586, 25704, 38085,
                                                 16944, 17304, 28584, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 40773, 0, 3, 35586, 36027, 26124, 38673,
                                                 17304, 17664, 29124, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41529, 3, 18384, 18399, 29685, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41557, 3, 18399, 18414, 29706, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41585, 3, 18414, 18429, 29727, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41613, 3, 18429, 18444, 29748, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41641, 3, 18444, 18459, 29769, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41669, 3, 18459, 18474, 29790, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41697, 3, 18474, 18489, 29811, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 41725, 0, 3, 29664, 41529, 29895, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41809, 0, 3, 29685, 41557, 29958, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41893, 0, 3, 29706, 41585, 30021, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41977, 0, 3, 29727, 41613, 30084, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 42061, 0, 3, 29748, 41641, 30147, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 42145, 0, 3, 29769, 41669, 30210, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 42229, 0, 3, 29790, 41697, 30273, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 42313, 0, 3, 29832, 41725, 18924, 19014,
                                                 30462, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42481, 0, 3, 29895, 41809, 19014, 19104,
                                                 30588, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42649, 0, 3, 29958, 41893, 19104, 19194,
                                                 30714, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42817, 0, 3, 30021, 41977, 19194, 19284,
                                                 30840, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42985, 0, 3, 30084, 42061, 19284, 19374,
                                                 30966, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 43153, 0, 3, 30147, 42145, 19374, 19464,
                                                 31092, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 43321, 0, 3, 30210, 42229, 19464, 19554,
                                                 31218, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 43489, 0, 3, 30336, 42313, 19734, 19884,
                                                 31344, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 43769, 0, 3, 30462, 42481, 19884, 20034,
                                                 31554, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44049, 0, 3, 30588, 42649, 20034, 20184,
                                                 31764, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44329, 0, 3, 30714, 42817, 20184, 20334,
                                                 31974, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44609, 0, 3, 30840, 42985, 20334, 20484,
                                                 32184, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44889, 0, 3, 30966, 43153, 20484, 20634,
                                                 32394, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 45169, 0, 3, 31092, 43321, 20634, 20784,
                                                 32604, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 45449, 0, 3, 42313, 42481, 31554, 44049,
                                                 21084, 21309, 33129, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 45869, 0, 3, 42481, 42649, 31764, 44329,
                                                 21309, 21534, 33444, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 46289, 0, 3, 42649, 42817, 31974, 44609,
                                                 21534, 21759, 33759, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 46709, 0, 3, 42817, 42985, 32184, 44889,
                                                 21759, 21984, 34074, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 47129, 0, 3, 42985, 43153, 32394, 45169,
                                                 21984, 22209, 34389, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 47549, 0, 3, 43489, 43769, 32814, 45449,
                                                 22659, 22974, 34704, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 48137, 0, 3, 43769, 44049, 33129, 45869,
                                                 22974, 23289, 35145, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 48725, 0, 3, 44049, 44329, 33444, 46289,
                                                 23289, 23604, 35586, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 49313, 0, 3, 44329, 44609, 33759, 46709,
                                                 23604, 23919, 36027, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 49901, 0, 3, 44609, 44889, 34074, 47129,
                                                 23919, 24234, 36468, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 50489, 0, 3, 45449, 45869, 35145, 48725,
                                                 24864, 25284, 37497, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 51273, 0, 3, 45869, 46289, 35586, 49313,
                                                 25284, 25704, 38085, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 52057, 0, 3, 46289, 46709, 36027, 49901,
                                                 25704, 26124, 38673, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 52841, 0, 3, 47549, 48137, 36909, 50489,
                                                 26964, 27504, 39261, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 53849, 0, 3, 48137, 48725, 37497, 51273,
                                                 27504, 28044, 40017, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 54857, 0, 3, 48725, 49313, 38085, 52057,
                                                 28044, 28584, 40773, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 55865, 3, 29664, 29685, 41557, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 55901, 3, 29685, 29706, 41585, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 55937, 3, 29706, 29727, 41613, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 55973, 3, 29727, 29748, 41641, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 56009, 3, 29748, 29769, 41669, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 56045, 3, 29769, 29790, 41697, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 56081, 0, 3, 41529, 55865, 41809, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 56189, 0, 3, 41557, 55901, 41893, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 56297, 0, 3, 41585, 55937, 41977, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 56405, 0, 3, 41613, 55973, 42061, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 56513, 0, 3, 41641, 56009, 42145, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 56621, 0, 3, 41669, 56045, 42229, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 56729, 0, 3, 41725, 56081, 30336, 30462,
                                                 42481, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 56945, 0, 3, 41809, 56189, 30462, 30588,
                                                 42649, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 57161, 0, 3, 41893, 56297, 30588, 30714,
                                                 42817, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 57377, 0, 3, 41977, 56405, 30714, 30840,
                                                 42985, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 57593, 0, 3, 42061, 56513, 30840, 30966,
                                                 43153, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 57809, 0, 3, 42145, 56621, 30966, 31092,
                                                 43321, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 58025, 0, 3, 42481, 56945, 31344, 31554,
                                                 44049, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 58385, 0, 3, 42649, 57161, 31554, 31764,
                                                 44329, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 58745, 0, 3, 42817, 57377, 31764, 31974,
                                                 44609, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 59105, 0, 3, 42985, 57593, 31974, 32184,
                                                 44889, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 59465, 0, 3, 43153, 57809, 32184, 32394,
                                                 45169, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 59825, 0, 3, 56729, 56945, 44049, 58385,
                                                 32814, 33129, 45869, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 60365, 0, 3, 56945, 57161, 44329, 58745,
                                                 33129, 33444, 46289, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 60905, 0, 3, 57161, 57377, 44609, 59105,
                                                 33444, 33759, 46709, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 61445, 0, 3, 57377, 57593, 44889, 59465,
                                                 33759, 34074, 47129, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 61985, 0, 3, 58025, 58385, 45869, 60365,
                                                 34704, 35145, 48725, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 62741, 0, 3, 58385, 58745, 46289, 60905,
                                                 35145, 35586, 49313, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 63497, 0, 3, 58745, 59105, 46709, 61445,
                                                 35586, 36027, 49901, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 64253, 0, 3, 59825, 60365, 48725, 62741,
                                                 36909, 37497, 51273, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 65261, 0, 3, 60365, 60905, 49313, 63497,
                                                 37497, 38085, 52057, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 66269, 0, 3, 61985, 62741, 51273, 65261,
                                                 39261, 40017, 54857, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 67565, 3, 41529, 41557, 55901, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 67610, 3, 41557, 41585, 55937, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 67655, 3, 41585, 41613, 55973, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 67700, 3, 41613, 41641, 56009, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 67745, 3, 41641, 41669, 56045, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 67790, 0, 3, 55865, 67565, 56189, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 67925, 0, 3, 55901, 67610, 56297, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 68060, 0, 3, 55937, 67655, 56405, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 68195, 0, 3, 55973, 67700, 56513, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 68330, 0, 3, 56009, 67745, 56621, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 68465, 0, 3, 56081, 67790, 42313, 42481,
                                                 56945, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 68735, 0, 3, 56189, 67925, 42481, 42649,
                                                 57161, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 69005, 0, 3, 56297, 68060, 42649, 42817,
                                                 57377, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 69275, 0, 3, 56405, 68195, 42817, 42985,
                                                 57593, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 69545, 0, 3, 56513, 68330, 42985, 43153,
                                                 57809, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 69815, 0, 3, 56729, 68465, 43489, 43769,
                                                 58025, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 70265, 0, 3, 56945, 68735, 43769, 44049,
                                                 58385, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 70715, 0, 3, 57161, 69005, 44049, 44329,
                                                 58745, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 71165, 0, 3, 57377, 69275, 44329, 44609,
                                                 59105, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 71615, 0, 3, 57593, 69545, 44609, 44889,
                                                 59465, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 72065, 0, 3, 68465, 68735, 58385, 70715,
                                                 45449, 45869, 60365, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 72740, 0, 3, 68735, 69005, 58745, 71165,
                                                 45869, 46289, 60905, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 73415, 0, 3, 69005, 69275, 59105, 71615,
                                                 46289, 46709, 61445, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 74090, 0, 3, 69815, 70265, 59825, 72065,
                                                 47549, 48137, 61985, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 75035, 0, 3, 70265, 70715, 60365, 72740,
                                                 48137, 48725, 62741, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 75980, 0, 3, 70715, 71165, 60905, 73415,
                                                 48725, 49313, 63497, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 76925, 0, 3, 72065, 72740, 62741, 75980,
                                                 50489, 51273, 65261, ncols, alpha, beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 78185, 0, 3, 74090, 75035, 64253, 76925,
                                                 52841, 53849, 66269, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 79805, 78185, 1620, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 81425, 79805, 36, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 81425, 17, nmax);
}

}  // namespace simdt2ceri
