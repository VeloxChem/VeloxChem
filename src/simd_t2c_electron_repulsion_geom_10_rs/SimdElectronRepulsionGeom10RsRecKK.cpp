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


#include "SimdElectronRepulsionGeom10RsRecKK.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
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
#include "SimdGeometryK1.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_kk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_kk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 177522, 169206, 7776, nvalues);

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
                                                9, 10, 11, 12, 13, 14, 15}, ncols, fj, mu,
                                                omega);

            simdfunc::compute_boys_function(buffer, coordinates, 22, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 80, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 83, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 86, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 89, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 92, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 95, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 98, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 101, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 104, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 107, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 110, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 113, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 116, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 119, 0, 37, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 7, 8, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 8, 9, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 9, 10, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 10, 11, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 11, 12, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 12, 13, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 13, 14, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 164, 0, 14, 15, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 170, 0, 15, 16, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 176, 0, 16, 17, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 182, 0, 17, 18, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 188, 0, 18, 19, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 194, 0, 19, 20, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 200, 0, 23, 24, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 206, 0, 24, 25, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 212, 0, 25, 26, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 218, 0, 26, 27, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 224, 0, 27, 28, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 230, 0, 28, 29, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 236, 0, 29, 30, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 242, 0, 30, 31, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 248, 0, 31, 32, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 254, 0, 32, 33, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 260, 0, 33, 34, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 266, 0, 34, 35, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 272, 0, 35, 36, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 38, 41, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 41, 44, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 44, 47, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 47, 50, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 50, 53, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 53, 56, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 56, 59, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 59, 62, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 62, 65, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 65, 68, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 68, 71, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 71, 74, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 398, 0, 80, 83, 206, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 408, 0, 83, 86, 212, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 418, 0, 86, 89, 218, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 428, 0, 89, 92, 224, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 438, 0, 92, 95, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 448, 0, 95, 98, 236, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 458, 0, 98, 101, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 468, 0, 101, 104, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 478, 0, 104, 107, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 488, 0, 107, 110, 260, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 498, 0, 110, 113, 266, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 508, 0, 113, 116, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 518, 0, 122, 128, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 533, 0, 128, 134, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 548, 0, 134, 140, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 563, 0, 140, 146, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 578, 0, 146, 152, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 593, 0, 152, 158, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 608, 0, 158, 164, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 623, 0, 164, 170, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 638, 0, 170, 176, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 653, 0, 176, 182, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 668, 0, 182, 188, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 683, 0, 200, 206, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 698, 0, 206, 212, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 713, 0, 212, 218, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 728, 0, 218, 224, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 743, 0, 224, 230, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 758, 0, 230, 236, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 773, 0, 236, 242, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 788, 0, 242, 248, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 803, 0, 248, 254, 488, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 818, 0, 254, 260, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 833, 0, 260, 266, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 848, 0, 278, 288, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 869, 0, 288, 298, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 890, 0, 298, 308, 563, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 911, 0, 308, 318, 578, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 932, 0, 318, 328, 593, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 953, 0, 328, 338, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 974, 0, 338, 348, 623, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 995, 0, 348, 358, 638, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1016, 0, 358, 368, 653, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1037, 0, 368, 378, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1058, 0, 398, 408, 698, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1079, 0, 408, 418, 713, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1100, 0, 418, 428, 728, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1121, 0, 428, 438, 743, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1142, 0, 438, 448, 758, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1163, 0, 448, 458, 773, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1184, 0, 458, 468, 788, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1205, 0, 468, 478, 803, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1226, 0, 478, 488, 818, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1247, 0, 488, 498, 833, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1268, 0, 518, 533, 869, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1296, 0, 533, 548, 890, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1324, 0, 548, 563, 911, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1352, 0, 563, 578, 932, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1380, 0, 578, 593, 953, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1408, 0, 593, 608, 974, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1436, 0, 608, 623, 995, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1464, 0, 623, 638, 1016, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1492, 0, 638, 653, 1037, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1520, 0, 683, 698, 1079, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1548, 0, 698, 713, 1100, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1576, 0, 713, 728, 1121, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1604, 0, 728, 743, 1142, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1632, 0, 743, 758, 1163, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1660, 0, 758, 773, 1184, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1688, 0, 773, 788, 1205, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1716, 0, 788, 803, 1226, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1744, 0, 803, 818, 1247, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1772, 0, 848, 869, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1808, 0, 869, 890, 1324, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1844, 0, 890, 911, 1352, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1880, 0, 911, 932, 1380, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1916, 0, 932, 953, 1408, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1952, 0, 953, 974, 1436, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1988, 0, 974, 995, 1464, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2024, 0, 995, 1016, 1492, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2060, 0, 1058, 1079, 1548, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2096, 0, 1079, 1100, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2132, 0, 1100, 1121, 1604, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2168, 0, 1121, 1142, 1632, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2204, 0, 1142, 1163, 1660, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2240, 0, 1163, 1184, 1688, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2276, 0, 1184, 1205, 1716, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2312, 0, 1205, 1226, 1744, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2348, 0, 1268, 1296, 1808, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2393, 0, 1296, 1324, 1844, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2438, 0, 1324, 1352, 1880, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2483, 0, 1352, 1380, 1916, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2528, 0, 1380, 1408, 1952, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2573, 0, 1408, 1436, 1988, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2618, 0, 1436, 1464, 2024, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2663, 0, 1520, 1548, 2096, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2708, 0, 1548, 1576, 2132, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2753, 0, 1576, 1604, 2168, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2798, 0, 1604, 1632, 2204, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2843, 0, 1632, 1660, 2240, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2888, 0, 1660, 1688, 2276, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2933, 0, 1688, 1716, 2312, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2978, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2981, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2984, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2987, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2990, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2993, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2996, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2999, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3002, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3005, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3008, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3011, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3014, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3017, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3020, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3023, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3026, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3029, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3032, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3035, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3038, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3041, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3044, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3047, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3050, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3053, 3, 37, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 3056, 3, 8, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3065, 3, 9, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3074, 3, 10, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3083, 3, 11, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3092, 3, 12, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3101, 3, 13, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3110, 3, 14, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3119, 3, 15, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3128, 3, 16, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3137, 3, 17, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3146, 3, 18, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3155, 3, 19, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3164, 3, 20, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3173, 3, 24, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3182, 3, 25, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3191, 3, 26, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3200, 3, 27, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3209, 3, 28, 95, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3218, 3, 29, 98, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3227, 3, 30, 101, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3236, 3, 31, 104, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3245, 3, 32, 107, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3254, 3, 33, 110, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3263, 3, 34, 113, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3272, 3, 35, 116, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3281, 3, 36, 119, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3290, 0, 3, 38, 3056, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3308, 0, 3, 41, 3065, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3326, 0, 3, 44, 3074, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3344, 0, 3, 47, 3083, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3362, 0, 3, 50, 3092, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3380, 0, 3, 53, 3101, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3398, 0, 3, 56, 3110, 158, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3416, 0, 3, 59, 3119, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3434, 0, 3, 62, 3128, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3452, 0, 3, 65, 3137, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3470, 0, 3, 68, 3146, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3488, 0, 3, 71, 3155, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3506, 0, 3, 74, 3164, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3524, 0, 3, 80, 3173, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3542, 0, 3, 83, 3182, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3560, 0, 3, 86, 3191, 212, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3578, 0, 3, 89, 3200, 218, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3596, 0, 3, 92, 3209, 224, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3614, 0, 3, 95, 3218, 230, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3632, 0, 3, 98, 3227, 236, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3650, 0, 3, 101, 3236, 242, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3668, 0, 3, 104, 3245, 248, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3686, 0, 3, 107, 3254, 254, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3704, 0, 3, 110, 3263, 260, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3722, 0, 3, 113, 3272, 266, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3740, 0, 3, 116, 3281, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3758, 0, 3, 128, 3326, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3788, 0, 3, 134, 3344, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3818, 0, 3, 140, 3362, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3848, 0, 3, 146, 3380, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3878, 0, 3, 152, 3398, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3908, 0, 3, 158, 3416, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3938, 0, 3, 164, 3434, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3968, 0, 3, 170, 3452, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3998, 0, 3, 176, 3470, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4028, 0, 3, 182, 3488, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4058, 0, 3, 188, 3506, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4088, 0, 3, 206, 3560, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4118, 0, 3, 212, 3578, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4148, 0, 3, 218, 3596, 428, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4178, 0, 3, 224, 3614, 438, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4208, 0, 3, 230, 3632, 448, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4238, 0, 3, 236, 3650, 458, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4268, 0, 3, 242, 3668, 468, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4298, 0, 3, 248, 3686, 478, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4328, 0, 3, 254, 3704, 488, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4358, 0, 3, 260, 3722, 498, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4388, 0, 3, 266, 3740, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4418, 0, 3, 278, 3758, 518, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4463, 0, 3, 288, 3788, 533, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4508, 0, 3, 298, 3818, 548, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4553, 0, 3, 308, 3848, 563, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4598, 0, 3, 318, 3878, 578, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4643, 0, 3, 328, 3908, 593, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4688, 0, 3, 338, 3938, 608, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4733, 0, 3, 348, 3968, 623, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4778, 0, 3, 358, 3998, 638, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4823, 0, 3, 368, 4028, 653, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4868, 0, 3, 378, 4058, 668, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4913, 0, 3, 398, 4088, 683, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4958, 0, 3, 408, 4118, 698, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5003, 0, 3, 418, 4148, 713, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5048, 0, 3, 428, 4178, 728, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5093, 0, 3, 438, 4208, 743, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5138, 0, 3, 448, 4238, 758, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5183, 0, 3, 458, 4268, 773, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5228, 0, 3, 468, 4298, 788, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5273, 0, 3, 478, 4328, 803, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5318, 0, 3, 488, 4358, 818, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5363, 0, 3, 498, 4388, 833, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5408, 0, 3, 533, 4508, 869, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5471, 0, 3, 548, 4553, 890, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5534, 0, 3, 563, 4598, 911, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5597, 0, 3, 578, 4643, 932, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5660, 0, 3, 593, 4688, 953, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5723, 0, 3, 608, 4733, 974, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5786, 0, 3, 623, 4778, 995, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5849, 0, 3, 638, 4823, 1016, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5912, 0, 3, 653, 4868, 1037, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5975, 0, 3, 698, 5003, 1079, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6038, 0, 3, 713, 5048, 1100, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6101, 0, 3, 728, 5093, 1121, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6164, 0, 3, 743, 5138, 1142, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6227, 0, 3, 758, 5183, 1163, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6290, 0, 3, 773, 5228, 1184, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6353, 0, 3, 788, 5273, 1205, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6416, 0, 3, 803, 5318, 1226, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6479, 0, 3, 818, 5363, 1247, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6542, 0, 3, 848, 5408, 1268, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6626, 0, 3, 869, 5471, 1296, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6710, 0, 3, 890, 5534, 1324, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6794, 0, 3, 911, 5597, 1352, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6878, 0, 3, 932, 5660, 1380, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6962, 0, 3, 953, 5723, 1408, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7046, 0, 3, 974, 5786, 1436, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7130, 0, 3, 995, 5849, 1464, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7214, 0, 3, 1016, 5912, 1492, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7298, 0, 3, 1058, 5975, 1520, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7382, 0, 3, 1079, 6038, 1548, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7466, 0, 3, 1100, 6101, 1576, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7550, 0, 3, 1121, 6164, 1604, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7634, 0, 3, 1142, 6227, 1632, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7718, 0, 3, 1163, 6290, 1660, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7802, 0, 3, 1184, 6353, 1688, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7886, 0, 3, 1205, 6416, 1716, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7970, 0, 3, 1226, 6479, 1744, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8054, 0, 3, 1296, 6710, 1808, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8162, 0, 3, 1324, 6794, 1844, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8270, 0, 3, 1352, 6878, 1880, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8378, 0, 3, 1380, 6962, 1916, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8486, 0, 3, 1408, 7046, 1952, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8594, 0, 3, 1436, 7130, 1988, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8702, 0, 3, 1464, 7214, 2024, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8810, 0, 3, 1548, 7466, 2096, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8918, 0, 3, 1576, 7550, 2132, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9026, 0, 3, 1604, 7634, 2168, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9134, 0, 3, 1632, 7718, 2204, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9242, 0, 3, 1660, 7802, 2240, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9350, 0, 3, 1688, 7886, 2276, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9458, 0, 3, 1716, 7970, 2312, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9566, 0, 3, 1772, 8054, 2348, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9701, 0, 3, 1808, 8162, 2393, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9836, 0, 3, 1844, 8270, 2438, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9971, 0, 3, 1880, 8378, 2483, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10106, 0, 3, 1916, 8486, 2528, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10241, 0, 3, 1952, 8594, 2573, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10376, 0, 3, 1988, 8702, 2618, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10511, 0, 3, 2060, 8810, 2663, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10646, 0, 3, 2096, 8918, 2708, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10781, 0, 3, 2132, 9026, 2753, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10916, 0, 3, 2168, 9134, 2798, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11051, 0, 3, 2204, 9242, 2843, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11186, 0, 3, 2240, 9350, 2888, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11321, 0, 3, 2276, 9458, 2933, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 11456, 3, 8, 9, 2981, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11462, 3, 9, 10, 2984, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11468, 3, 10, 11, 2987, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11474, 3, 11, 12, 2990, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11480, 3, 12, 13, 2993, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11486, 3, 13, 14, 2996, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11492, 3, 14, 15, 2999, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11498, 3, 15, 16, 3002, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11504, 3, 16, 17, 3005, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11510, 3, 17, 18, 3008, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11516, 3, 18, 19, 3011, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11522, 3, 19, 20, 3014, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11528, 3, 24, 25, 3020, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11534, 3, 25, 26, 3023, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11540, 3, 26, 27, 3026, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11546, 3, 27, 28, 3029, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11552, 3, 28, 29, 3032, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11558, 3, 29, 30, 3035, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11564, 3, 30, 31, 3038, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11570, 3, 31, 32, 3041, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11576, 3, 32, 33, 3044, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11582, 3, 33, 34, 3047, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11588, 3, 34, 35, 3050, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 11594, 3, 35, 36, 3053, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 11600, 0, 3, 2978, 11456, 3065, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11618, 0, 3, 2981, 11462, 3074, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11636, 0, 3, 2984, 11468, 3083, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11654, 0, 3, 2987, 11474, 3092, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11672, 0, 3, 2990, 11480, 3101, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11690, 0, 3, 2993, 11486, 3110, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11708, 0, 3, 2996, 11492, 3119, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11726, 0, 3, 2999, 11498, 3128, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11744, 0, 3, 3002, 11504, 3137, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11762, 0, 3, 3005, 11510, 3146, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11780, 0, 3, 3008, 11516, 3155, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11798, 0, 3, 3011, 11522, 3164, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11816, 0, 3, 3017, 11528, 3182, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11834, 0, 3, 3020, 11534, 3191, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11852, 0, 3, 3023, 11540, 3200, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11870, 0, 3, 3026, 11546, 3209, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11888, 0, 3, 3029, 11552, 3218, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11906, 0, 3, 3032, 11558, 3227, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11924, 0, 3, 3035, 11564, 3236, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11942, 0, 3, 3038, 11570, 3245, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11960, 0, 3, 3041, 11576, 3254, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11978, 0, 3, 3044, 11582, 3263, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 11996, 0, 3, 3047, 11588, 3272, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12014, 0, 3, 3050, 11594, 3281, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 12032, 0, 3, 3065, 11618, 122, 128,
                                                 3326, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12068, 0, 3, 3074, 11636, 128, 134,
                                                 3344, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12104, 0, 3, 3083, 11654, 134, 140,
                                                 3362, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12140, 0, 3, 3092, 11672, 140, 146,
                                                 3380, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12176, 0, 3, 3101, 11690, 146, 152,
                                                 3398, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12212, 0, 3, 3110, 11708, 152, 158,
                                                 3416, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12248, 0, 3, 3119, 11726, 158, 164,
                                                 3434, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12284, 0, 3, 3128, 11744, 164, 170,
                                                 3452, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12320, 0, 3, 3137, 11762, 170, 176,
                                                 3470, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12356, 0, 3, 3146, 11780, 176, 182,
                                                 3488, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12392, 0, 3, 3155, 11798, 182, 188,
                                                 3506, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12428, 0, 3, 3182, 11834, 200, 206,
                                                 3560, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12464, 0, 3, 3191, 11852, 206, 212,
                                                 3578, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12500, 0, 3, 3200, 11870, 212, 218,
                                                 3596, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12536, 0, 3, 3209, 11888, 218, 224,
                                                 3614, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12572, 0, 3, 3218, 11906, 224, 230,
                                                 3632, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12608, 0, 3, 3227, 11924, 230, 236,
                                                 3650, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12644, 0, 3, 3236, 11942, 236, 242,
                                                 3668, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12680, 0, 3, 3245, 11960, 242, 248,
                                                 3686, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12716, 0, 3, 3254, 11978, 248, 254,
                                                 3704, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12752, 0, 3, 3263, 11996, 254, 260,
                                                 3722, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12788, 0, 3, 3272, 12014, 260, 266,
                                                 3740, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12824, 0, 3, 3326, 12068, 278, 288,
                                                 3788, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12884, 0, 3, 3344, 12104, 288, 298,
                                                 3818, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12944, 0, 3, 3362, 12140, 298, 308,
                                                 3848, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13004, 0, 3, 3380, 12176, 308, 318,
                                                 3878, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13064, 0, 3, 3398, 12212, 318, 328,
                                                 3908, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13124, 0, 3, 3416, 12248, 328, 338,
                                                 3938, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13184, 0, 3, 3434, 12284, 338, 348,
                                                 3968, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13244, 0, 3, 3452, 12320, 348, 358,
                                                 3998, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13304, 0, 3, 3470, 12356, 358, 368,
                                                 4028, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13364, 0, 3, 3488, 12392, 368, 378,
                                                 4058, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13424, 0, 3, 3560, 12464, 398, 408,
                                                 4118, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13484, 0, 3, 3578, 12500, 408, 418,
                                                 4148, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13544, 0, 3, 3596, 12536, 418, 428,
                                                 4178, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13604, 0, 3, 3614, 12572, 428, 438,
                                                 4208, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13664, 0, 3, 3632, 12608, 438, 448,
                                                 4238, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13724, 0, 3, 3650, 12644, 448, 458,
                                                 4268, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13784, 0, 3, 3668, 12680, 458, 468,
                                                 4298, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13844, 0, 3, 3686, 12716, 468, 478,
                                                 4328, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13904, 0, 3, 3704, 12752, 478, 488,
                                                 4358, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13964, 0, 3, 3722, 12788, 488, 498,
                                                 4388, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14024, 0, 3, 12032, 12068, 3788, 12884,
                                                 518, 533, 4508, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14114, 0, 3, 12068, 12104, 3818, 12944,
                                                 533, 548, 4553, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14204, 0, 3, 12104, 12140, 3848, 13004,
                                                 548, 563, 4598, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14294, 0, 3, 12140, 12176, 3878, 13064,
                                                 563, 578, 4643, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14384, 0, 3, 12176, 12212, 3908, 13124,
                                                 578, 593, 4688, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14474, 0, 3, 12212, 12248, 3938, 13184,
                                                 593, 608, 4733, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14564, 0, 3, 12248, 12284, 3968, 13244,
                                                 608, 623, 4778, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14654, 0, 3, 12284, 12320, 3998, 13304,
                                                 623, 638, 4823, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14744, 0, 3, 12320, 12356, 4028, 13364,
                                                 638, 653, 4868, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14834, 0, 3, 12428, 12464, 4118, 13484,
                                                 683, 698, 5003, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14924, 0, 3, 12464, 12500, 4148, 13544,
                                                 698, 713, 5048, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15014, 0, 3, 12500, 12536, 4178, 13604,
                                                 713, 728, 5093, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15104, 0, 3, 12536, 12572, 4208, 13664,
                                                 728, 743, 5138, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15194, 0, 3, 12572, 12608, 4238, 13724,
                                                 743, 758, 5183, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15284, 0, 3, 12608, 12644, 4268, 13784,
                                                 758, 773, 5228, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15374, 0, 3, 12644, 12680, 4298, 13844,
                                                 773, 788, 5273, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15464, 0, 3, 12680, 12716, 4328, 13904,
                                                 788, 803, 5318, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15554, 0, 3, 12716, 12752, 4358, 13964,
                                                 803, 818, 5363, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15644, 0, 3, 12824, 12884, 4508, 14114,
                                                 848, 869, 5471, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15770, 0, 3, 12884, 12944, 4553, 14204,
                                                 869, 890, 5534, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15896, 0, 3, 12944, 13004, 4598, 14294,
                                                 890, 911, 5597, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16022, 0, 3, 13004, 13064, 4643, 14384,
                                                 911, 932, 5660, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16148, 0, 3, 13064, 13124, 4688, 14474,
                                                 932, 953, 5723, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16274, 0, 3, 13124, 13184, 4733, 14564,
                                                 953, 974, 5786, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16400, 0, 3, 13184, 13244, 4778, 14654,
                                                 974, 995, 5849, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16526, 0, 3, 13244, 13304, 4823, 14744,
                                                 995, 1016, 5912, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16652, 0, 3, 13424, 13484, 5003, 14924,
                                                 1058, 1079, 6038, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16778, 0, 3, 13484, 13544, 5048, 15014,
                                                 1079, 1100, 6101, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16904, 0, 3, 13544, 13604, 5093, 15104,
                                                 1100, 1121, 6164, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17030, 0, 3, 13604, 13664, 5138, 15194,
                                                 1121, 1142, 6227, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17156, 0, 3, 13664, 13724, 5183, 15284,
                                                 1142, 1163, 6290, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17282, 0, 3, 13724, 13784, 5228, 15374,
                                                 1163, 1184, 6353, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17408, 0, 3, 13784, 13844, 5273, 15464,
                                                 1184, 1205, 6416, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17534, 0, 3, 13844, 13904, 5318, 15554,
                                                 1205, 1226, 6479, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17660, 0, 3, 14024, 14114, 5471, 15770,
                                                 1268, 1296, 6710, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17828, 0, 3, 14114, 14204, 5534, 15896,
                                                 1296, 1324, 6794, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17996, 0, 3, 14204, 14294, 5597, 16022,
                                                 1324, 1352, 6878, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18164, 0, 3, 14294, 14384, 5660, 16148,
                                                 1352, 1380, 6962, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18332, 0, 3, 14384, 14474, 5723, 16274,
                                                 1380, 1408, 7046, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18500, 0, 3, 14474, 14564, 5786, 16400,
                                                 1408, 1436, 7130, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18668, 0, 3, 14564, 14654, 5849, 16526,
                                                 1436, 1464, 7214, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18836, 0, 3, 14834, 14924, 6038, 16778,
                                                 1520, 1548, 7466, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19004, 0, 3, 14924, 15014, 6101, 16904,
                                                 1548, 1576, 7550, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19172, 0, 3, 15014, 15104, 6164, 17030,
                                                 1576, 1604, 7634, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19340, 0, 3, 15104, 15194, 6227, 17156,
                                                 1604, 1632, 7718, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19508, 0, 3, 15194, 15284, 6290, 17282,
                                                 1632, 1660, 7802, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19676, 0, 3, 15284, 15374, 6353, 17408,
                                                 1660, 1688, 7886, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19844, 0, 3, 15374, 15464, 6416, 17534,
                                                 1688, 1716, 7970, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20012, 0, 3, 15644, 15770, 6710, 17828,
                                                 1772, 1808, 8162, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20228, 0, 3, 15770, 15896, 6794, 17996,
                                                 1808, 1844, 8270, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20444, 0, 3, 15896, 16022, 6878, 18164,
                                                 1844, 1880, 8378, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20660, 0, 3, 16022, 16148, 6962, 18332,
                                                 1880, 1916, 8486, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20876, 0, 3, 16148, 16274, 7046, 18500,
                                                 1916, 1952, 8594, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21092, 0, 3, 16274, 16400, 7130, 18668,
                                                 1952, 1988, 8702, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21308, 0, 3, 16652, 16778, 7466, 19004,
                                                 2060, 2096, 8918, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21524, 0, 3, 16778, 16904, 7550, 19172,
                                                 2096, 2132, 9026, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21740, 0, 3, 16904, 17030, 7634, 19340,
                                                 2132, 2168, 9134, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21956, 0, 3, 17030, 17156, 7718, 19508,
                                                 2168, 2204, 9242, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22172, 0, 3, 17156, 17282, 7802, 19676,
                                                 2204, 2240, 9350, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22388, 0, 3, 17282, 17408, 7886, 19844,
                                                 2240, 2276, 9458, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22604, 0, 3, 17660, 17828, 8162, 20228,
                                                 2348, 2393, 9836, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22874, 0, 3, 17828, 17996, 8270, 20444,
                                                 2393, 2438, 9971, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23144, 0, 3, 17996, 18164, 8378, 20660,
                                                 2438, 2483, 10106, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23414, 0, 3, 18164, 18332, 8486, 20876,
                                                 2483, 2528, 10241, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23684, 0, 3, 18332, 18500, 8594, 21092,
                                                 2528, 2573, 10376, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23954, 0, 3, 18836, 19004, 8918, 21524,
                                                 2663, 2708, 10781, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 24224, 0, 3, 19004, 19172, 9026, 21740,
                                                 2708, 2753, 10916, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 24494, 0, 3, 19172, 19340, 9134, 21956,
                                                 2753, 2798, 11051, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 24764, 0, 3, 19340, 19508, 9242, 22172,
                                                 2798, 2843, 11186, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25034, 0, 3, 19508, 19676, 9350, 22388,
                                                 2843, 2888, 11321, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25304, 3, 2978, 2981, 11462, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25314, 3, 2981, 2984, 11468, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25324, 3, 2984, 2987, 11474, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25334, 3, 2987, 2990, 11480, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25344, 3, 2990, 2993, 11486, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25354, 3, 2993, 2996, 11492, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25364, 3, 2996, 2999, 11498, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25374, 3, 2999, 3002, 11504, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25384, 3, 3002, 3005, 11510, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25394, 3, 3005, 3008, 11516, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25404, 3, 3008, 3011, 11522, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25414, 3, 3017, 3020, 11534, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25424, 3, 3020, 3023, 11540, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25434, 3, 3023, 3026, 11546, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25444, 3, 3026, 3029, 11552, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25454, 3, 3029, 3032, 11558, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25464, 3, 3032, 3035, 11564, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25474, 3, 3035, 3038, 11570, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25484, 3, 3038, 3041, 11576, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25494, 3, 3041, 3044, 11582, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25504, 3, 3044, 3047, 11588, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25514, 3, 3047, 3050, 11594, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 25524, 0, 3, 11456, 25304, 11618, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25554, 0, 3, 11462, 25314, 11636, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25584, 0, 3, 11468, 25324, 11654, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25614, 0, 3, 11474, 25334, 11672, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25644, 0, 3, 11480, 25344, 11690, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25674, 0, 3, 11486, 25354, 11708, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25704, 0, 3, 11492, 25364, 11726, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25734, 0, 3, 11498, 25374, 11744, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25764, 0, 3, 11504, 25384, 11762, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25794, 0, 3, 11510, 25394, 11780, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25824, 0, 3, 11516, 25404, 11798, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25854, 0, 3, 11528, 25414, 11834, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25884, 0, 3, 11534, 25424, 11852, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25914, 0, 3, 11540, 25434, 11870, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25944, 0, 3, 11546, 25444, 11888, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 25974, 0, 3, 11552, 25454, 11906, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26004, 0, 3, 11558, 25464, 11924, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26034, 0, 3, 11564, 25474, 11942, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26064, 0, 3, 11570, 25484, 11960, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26094, 0, 3, 11576, 25494, 11978, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26124, 0, 3, 11582, 25504, 11996, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26154, 0, 3, 11588, 25514, 12014, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 26184, 0, 3, 11600, 25524, 3290, 3308,
                                                 12032, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26244, 0, 3, 11618, 25554, 3308, 3326,
                                                 12068, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26304, 0, 3, 11636, 25584, 3326, 3344,
                                                 12104, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26364, 0, 3, 11654, 25614, 3344, 3362,
                                                 12140, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26424, 0, 3, 11672, 25644, 3362, 3380,
                                                 12176, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26484, 0, 3, 11690, 25674, 3380, 3398,
                                                 12212, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26544, 0, 3, 11708, 25704, 3398, 3416,
                                                 12248, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26604, 0, 3, 11726, 25734, 3416, 3434,
                                                 12284, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26664, 0, 3, 11744, 25764, 3434, 3452,
                                                 12320, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26724, 0, 3, 11762, 25794, 3452, 3470,
                                                 12356, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26784, 0, 3, 11780, 25824, 3470, 3488,
                                                 12392, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26844, 0, 3, 11816, 25854, 3524, 3542,
                                                 12428, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26904, 0, 3, 11834, 25884, 3542, 3560,
                                                 12464, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26964, 0, 3, 11852, 25914, 3560, 3578,
                                                 12500, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27024, 0, 3, 11870, 25944, 3578, 3596,
                                                 12536, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27084, 0, 3, 11888, 25974, 3596, 3614,
                                                 12572, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27144, 0, 3, 11906, 26004, 3614, 3632,
                                                 12608, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27204, 0, 3, 11924, 26034, 3632, 3650,
                                                 12644, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27264, 0, 3, 11942, 26064, 3650, 3668,
                                                 12680, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27324, 0, 3, 11960, 26094, 3668, 3686,
                                                 12716, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27384, 0, 3, 11978, 26124, 3686, 3704,
                                                 12752, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27444, 0, 3, 11996, 26154, 3704, 3722,
                                                 12788, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27504, 0, 3, 12068, 26304, 3758, 3788,
                                                 12884, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27604, 0, 3, 12104, 26364, 3788, 3818,
                                                 12944, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27704, 0, 3, 12140, 26424, 3818, 3848,
                                                 13004, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27804, 0, 3, 12176, 26484, 3848, 3878,
                                                 13064, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27904, 0, 3, 12212, 26544, 3878, 3908,
                                                 13124, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28004, 0, 3, 12248, 26604, 3908, 3938,
                                                 13184, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28104, 0, 3, 12284, 26664, 3938, 3968,
                                                 13244, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28204, 0, 3, 12320, 26724, 3968, 3998,
                                                 13304, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28304, 0, 3, 12356, 26784, 3998, 4028,
                                                 13364, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28404, 0, 3, 12464, 26964, 4088, 4118,
                                                 13484, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28504, 0, 3, 12500, 27024, 4118, 4148,
                                                 13544, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28604, 0, 3, 12536, 27084, 4148, 4178,
                                                 13604, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28704, 0, 3, 12572, 27144, 4178, 4208,
                                                 13664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28804, 0, 3, 12608, 27204, 4208, 4238,
                                                 13724, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28904, 0, 3, 12644, 27264, 4238, 4268,
                                                 13784, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29004, 0, 3, 12680, 27324, 4268, 4298,
                                                 13844, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29104, 0, 3, 12716, 27384, 4298, 4328,
                                                 13904, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29204, 0, 3, 12752, 27444, 4328, 4358,
                                                 13964, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29304, 0, 3, 26184, 26244, 12824, 27504,
                                                 4418, 4463, 14024, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29454, 0, 3, 26244, 26304, 12884, 27604,
                                                 4463, 4508, 14114, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29604, 0, 3, 26304, 26364, 12944, 27704,
                                                 4508, 4553, 14204, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29754, 0, 3, 26364, 26424, 13004, 27804,
                                                 4553, 4598, 14294, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29904, 0, 3, 26424, 26484, 13064, 27904,
                                                 4598, 4643, 14384, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30054, 0, 3, 26484, 26544, 13124, 28004,
                                                 4643, 4688, 14474, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30204, 0, 3, 26544, 26604, 13184, 28104,
                                                 4688, 4733, 14564, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30354, 0, 3, 26604, 26664, 13244, 28204,
                                                 4733, 4778, 14654, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30504, 0, 3, 26664, 26724, 13304, 28304,
                                                 4778, 4823, 14744, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30654, 0, 3, 26844, 26904, 13424, 28404,
                                                 4913, 4958, 14834, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30804, 0, 3, 26904, 26964, 13484, 28504,
                                                 4958, 5003, 14924, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30954, 0, 3, 26964, 27024, 13544, 28604,
                                                 5003, 5048, 15014, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31104, 0, 3, 27024, 27084, 13604, 28704,
                                                 5048, 5093, 15104, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31254, 0, 3, 27084, 27144, 13664, 28804,
                                                 5093, 5138, 15194, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31404, 0, 3, 27144, 27204, 13724, 28904,
                                                 5138, 5183, 15284, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31554, 0, 3, 27204, 27264, 13784, 29004,
                                                 5183, 5228, 15374, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31704, 0, 3, 27264, 27324, 13844, 29104,
                                                 5228, 5273, 15464, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31854, 0, 3, 27324, 27384, 13904, 29204,
                                                 5273, 5318, 15554, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32004, 0, 3, 27504, 27604, 14114, 29604,
                                                 5408, 5471, 15770, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32214, 0, 3, 27604, 27704, 14204, 29754,
                                                 5471, 5534, 15896, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32424, 0, 3, 27704, 27804, 14294, 29904,
                                                 5534, 5597, 16022, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32634, 0, 3, 27804, 27904, 14384, 30054,
                                                 5597, 5660, 16148, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32844, 0, 3, 27904, 28004, 14474, 30204,
                                                 5660, 5723, 16274, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33054, 0, 3, 28004, 28104, 14564, 30354,
                                                 5723, 5786, 16400, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33264, 0, 3, 28104, 28204, 14654, 30504,
                                                 5786, 5849, 16526, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33474, 0, 3, 28404, 28504, 14924, 30954,
                                                 5975, 6038, 16778, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33684, 0, 3, 28504, 28604, 15014, 31104,
                                                 6038, 6101, 16904, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33894, 0, 3, 28604, 28704, 15104, 31254,
                                                 6101, 6164, 17030, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34104, 0, 3, 28704, 28804, 15194, 31404,
                                                 6164, 6227, 17156, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34314, 0, 3, 28804, 28904, 15284, 31554,
                                                 6227, 6290, 17282, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34524, 0, 3, 28904, 29004, 15374, 31704,
                                                 6290, 6353, 17408, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34734, 0, 3, 29004, 29104, 15464, 31854,
                                                 6353, 6416, 17534, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 34944, 0, 3, 29304, 29454, 15644, 32004,
                                                 6542, 6626, 17660, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35224, 0, 3, 29454, 29604, 15770, 32214,
                                                 6626, 6710, 17828, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35504, 0, 3, 29604, 29754, 15896, 32424,
                                                 6710, 6794, 17996, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35784, 0, 3, 29754, 29904, 16022, 32634,
                                                 6794, 6878, 18164, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36064, 0, 3, 29904, 30054, 16148, 32844,
                                                 6878, 6962, 18332, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36344, 0, 3, 30054, 30204, 16274, 33054,
                                                 6962, 7046, 18500, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36624, 0, 3, 30204, 30354, 16400, 33264,
                                                 7046, 7130, 18668, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36904, 0, 3, 30654, 30804, 16652, 33474,
                                                 7298, 7382, 18836, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 37184, 0, 3, 30804, 30954, 16778, 33684,
                                                 7382, 7466, 19004, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 37464, 0, 3, 30954, 31104, 16904, 33894,
                                                 7466, 7550, 19172, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 37744, 0, 3, 31104, 31254, 17030, 34104,
                                                 7550, 7634, 19340, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 38024, 0, 3, 31254, 31404, 17156, 34314,
                                                 7634, 7718, 19508, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 38304, 0, 3, 31404, 31554, 17282, 34524,
                                                 7718, 7802, 19676, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 38584, 0, 3, 31554, 31704, 17408, 34734,
                                                 7802, 7886, 19844, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 38864, 0, 3, 32004, 32214, 17828, 35504,
                                                 8054, 8162, 20228, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39224, 0, 3, 32214, 32424, 17996, 35784,
                                                 8162, 8270, 20444, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39584, 0, 3, 32424, 32634, 18164, 36064,
                                                 8270, 8378, 20660, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39944, 0, 3, 32634, 32844, 18332, 36344,
                                                 8378, 8486, 20876, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 40304, 0, 3, 32844, 33054, 18500, 36624,
                                                 8486, 8594, 21092, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 40664, 0, 3, 33474, 33684, 19004, 37464,
                                                 8810, 8918, 21524, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 41024, 0, 3, 33684, 33894, 19172, 37744,
                                                 8918, 9026, 21740, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 41384, 0, 3, 33894, 34104, 19340, 38024,
                                                 9026, 9134, 21956, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 41744, 0, 3, 34104, 34314, 19508, 38304,
                                                 9134, 9242, 22172, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 42104, 0, 3, 34314, 34524, 19676, 38584,
                                                 9242, 9350, 22388, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 42464, 0, 3, 34944, 35224, 20012, 38864,
                                                 9566, 9701, 22604, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 42914, 0, 3, 35224, 35504, 20228, 39224,
                                                 9701, 9836, 22874, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 43364, 0, 3, 35504, 35784, 20444, 39584,
                                                 9836, 9971, 23144, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 43814, 0, 3, 35784, 36064, 20660, 39944,
                                                 9971, 10106, 23414, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 44264, 0, 3, 36064, 36344, 20876, 40304,
                                                 10106, 10241, 23684, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 44714, 0, 3, 36904, 37184, 21308, 40664,
                                                 10511, 10646, 23954, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 45164, 0, 3, 37184, 37464, 21524, 41024,
                                                 10646, 10781, 24224, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 45614, 0, 3, 37464, 37744, 21740, 41384,
                                                 10781, 10916, 24494, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 46064, 0, 3, 37744, 38024, 21956, 41744,
                                                 10916, 11051, 24764, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 46514, 0, 3, 38024, 38304, 22172, 42104,
                                                 11051, 11186, 25034, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 46964, 3, 11456, 11462, 25314, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 46979, 3, 11462, 11468, 25324, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 46994, 3, 11468, 11474, 25334, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47009, 3, 11474, 11480, 25344, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47024, 3, 11480, 11486, 25354, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47039, 3, 11486, 11492, 25364, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47054, 3, 11492, 11498, 25374, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47069, 3, 11498, 11504, 25384, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47084, 3, 11504, 11510, 25394, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47099, 3, 11510, 11516, 25404, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47114, 3, 11528, 11534, 25424, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47129, 3, 11534, 11540, 25434, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47144, 3, 11540, 11546, 25444, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47159, 3, 11546, 11552, 25454, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47174, 3, 11552, 11558, 25464, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47189, 3, 11558, 11564, 25474, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47204, 3, 11564, 11570, 25484, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47219, 3, 11570, 11576, 25494, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47234, 3, 11576, 11582, 25504, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47249, 3, 11582, 11588, 25514, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 47264, 0, 3, 25304, 46964, 25554, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47309, 0, 3, 25314, 46979, 25584, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47354, 0, 3, 25324, 46994, 25614, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47399, 0, 3, 25334, 47009, 25644, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47444, 0, 3, 25344, 47024, 25674, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47489, 0, 3, 25354, 47039, 25704, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47534, 0, 3, 25364, 47054, 25734, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47579, 0, 3, 25374, 47069, 25764, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47624, 0, 3, 25384, 47084, 25794, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47669, 0, 3, 25394, 47099, 25824, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47714, 0, 3, 25414, 47114, 25884, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47759, 0, 3, 25424, 47129, 25914, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47804, 0, 3, 25434, 47144, 25944, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47849, 0, 3, 25444, 47159, 25974, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47894, 0, 3, 25454, 47174, 26004, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47939, 0, 3, 25464, 47189, 26034, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47984, 0, 3, 25474, 47204, 26064, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 48029, 0, 3, 25484, 47219, 26094, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 48074, 0, 3, 25494, 47234, 26124, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 48119, 0, 3, 25504, 47249, 26154, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 48164, 0, 3, 25554, 47309, 12032, 12068,
                                                 26304, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48254, 0, 3, 25584, 47354, 12068, 12104,
                                                 26364, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48344, 0, 3, 25614, 47399, 12104, 12140,
                                                 26424, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48434, 0, 3, 25644, 47444, 12140, 12176,
                                                 26484, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48524, 0, 3, 25674, 47489, 12176, 12212,
                                                 26544, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48614, 0, 3, 25704, 47534, 12212, 12248,
                                                 26604, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48704, 0, 3, 25734, 47579, 12248, 12284,
                                                 26664, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48794, 0, 3, 25764, 47624, 12284, 12320,
                                                 26724, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48884, 0, 3, 25794, 47669, 12320, 12356,
                                                 26784, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48974, 0, 3, 25884, 47759, 12428, 12464,
                                                 26964, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49064, 0, 3, 25914, 47804, 12464, 12500,
                                                 27024, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49154, 0, 3, 25944, 47849, 12500, 12536,
                                                 27084, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49244, 0, 3, 25974, 47894, 12536, 12572,
                                                 27144, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49334, 0, 3, 26004, 47939, 12572, 12608,
                                                 27204, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49424, 0, 3, 26034, 47984, 12608, 12644,
                                                 27264, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49514, 0, 3, 26064, 48029, 12644, 12680,
                                                 27324, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49604, 0, 3, 26094, 48074, 12680, 12716,
                                                 27384, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49694, 0, 3, 26124, 48119, 12716, 12752,
                                                 27444, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49784, 0, 3, 26304, 48254, 12824, 12884,
                                                 27604, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49934, 0, 3, 26364, 48344, 12884, 12944,
                                                 27704, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50084, 0, 3, 26424, 48434, 12944, 13004,
                                                 27804, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50234, 0, 3, 26484, 48524, 13004, 13064,
                                                 27904, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50384, 0, 3, 26544, 48614, 13064, 13124,
                                                 28004, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50534, 0, 3, 26604, 48704, 13124, 13184,
                                                 28104, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50684, 0, 3, 26664, 48794, 13184, 13244,
                                                 28204, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50834, 0, 3, 26724, 48884, 13244, 13304,
                                                 28304, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50984, 0, 3, 26964, 49064, 13424, 13484,
                                                 28504, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51134, 0, 3, 27024, 49154, 13484, 13544,
                                                 28604, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51284, 0, 3, 27084, 49244, 13544, 13604,
                                                 28704, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51434, 0, 3, 27144, 49334, 13604, 13664,
                                                 28804, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51584, 0, 3, 27204, 49424, 13664, 13724,
                                                 28904, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51734, 0, 3, 27264, 49514, 13724, 13784,
                                                 29004, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51884, 0, 3, 27324, 49604, 13784, 13844,
                                                 29104, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 52034, 0, 3, 27384, 49694, 13844, 13904,
                                                 29204, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52184, 0, 3, 48164, 48254, 27604, 49934,
                                                 14024, 14114, 29604, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52409, 0, 3, 48254, 48344, 27704, 50084,
                                                 14114, 14204, 29754, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52634, 0, 3, 48344, 48434, 27804, 50234,
                                                 14204, 14294, 29904, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52859, 0, 3, 48434, 48524, 27904, 50384,
                                                 14294, 14384, 30054, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53084, 0, 3, 48524, 48614, 28004, 50534,
                                                 14384, 14474, 30204, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53309, 0, 3, 48614, 48704, 28104, 50684,
                                                 14474, 14564, 30354, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53534, 0, 3, 48704, 48794, 28204, 50834,
                                                 14564, 14654, 30504, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53759, 0, 3, 48974, 49064, 28504, 51134,
                                                 14834, 14924, 30954, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53984, 0, 3, 49064, 49154, 28604, 51284,
                                                 14924, 15014, 31104, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 54209, 0, 3, 49154, 49244, 28704, 51434,
                                                 15014, 15104, 31254, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 54434, 0, 3, 49244, 49334, 28804, 51584,
                                                 15104, 15194, 31404, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 54659, 0, 3, 49334, 49424, 28904, 51734,
                                                 15194, 15284, 31554, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 54884, 0, 3, 49424, 49514, 29004, 51884,
                                                 15284, 15374, 31704, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 55109, 0, 3, 49514, 49604, 29104, 52034,
                                                 15374, 15464, 31854, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55334, 0, 3, 49784, 49934, 29604, 52409,
                                                 15644, 15770, 32214, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55649, 0, 3, 49934, 50084, 29754, 52634,
                                                 15770, 15896, 32424, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55964, 0, 3, 50084, 50234, 29904, 52859,
                                                 15896, 16022, 32634, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56279, 0, 3, 50234, 50384, 30054, 53084,
                                                 16022, 16148, 32844, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56594, 0, 3, 50384, 50534, 30204, 53309,
                                                 16148, 16274, 33054, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56909, 0, 3, 50534, 50684, 30354, 53534,
                                                 16274, 16400, 33264, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 57224, 0, 3, 50984, 51134, 30954, 53984,
                                                 16652, 16778, 33684, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 57539, 0, 3, 51134, 51284, 31104, 54209,
                                                 16778, 16904, 33894, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 57854, 0, 3, 51284, 51434, 31254, 54434,
                                                 16904, 17030, 34104, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 58169, 0, 3, 51434, 51584, 31404, 54659,
                                                 17030, 17156, 34314, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 58484, 0, 3, 51584, 51734, 31554, 54884,
                                                 17156, 17282, 34524, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 58799, 0, 3, 51734, 51884, 31704, 55109,
                                                 17282, 17408, 34734, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59114, 0, 3, 52184, 52409, 32214, 55649,
                                                 17660, 17828, 35504, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59534, 0, 3, 52409, 52634, 32424, 55964,
                                                 17828, 17996, 35784, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59954, 0, 3, 52634, 52859, 32634, 56279,
                                                 17996, 18164, 36064, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 60374, 0, 3, 52859, 53084, 32844, 56594,
                                                 18164, 18332, 36344, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 60794, 0, 3, 53084, 53309, 33054, 56909,
                                                 18332, 18500, 36624, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 61214, 0, 3, 53759, 53984, 33684, 57539,
                                                 18836, 19004, 37464, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 61634, 0, 3, 53984, 54209, 33894, 57854,
                                                 19004, 19172, 37744, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 62054, 0, 3, 54209, 54434, 34104, 58169,
                                                 19172, 19340, 38024, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 62474, 0, 3, 54434, 54659, 34314, 58484,
                                                 19340, 19508, 38304, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 62894, 0, 3, 54659, 54884, 34524, 58799,
                                                 19508, 19676, 38584, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 63314, 0, 3, 55334, 55649, 35504, 59534,
                                                 20012, 20228, 39224, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 63854, 0, 3, 55649, 55964, 35784, 59954,
                                                 20228, 20444, 39584, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 64394, 0, 3, 55964, 56279, 36064, 60374,
                                                 20444, 20660, 39944, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 64934, 0, 3, 56279, 56594, 36344, 60794,
                                                 20660, 20876, 40304, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 65474, 0, 3, 57224, 57539, 37464, 61634,
                                                 21308, 21524, 41024, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 66014, 0, 3, 57539, 57854, 37744, 62054,
                                                 21524, 21740, 41384, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 66554, 0, 3, 57854, 58169, 38024, 62474,
                                                 21740, 21956, 41744, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 67094, 0, 3, 58169, 58484, 38304, 62894,
                                                 21956, 22172, 42104, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 67634, 0, 3, 59114, 59534, 39224, 63854,
                                                 22604, 22874, 43364, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 68309, 0, 3, 59534, 59954, 39584, 64394,
                                                 22874, 23144, 43814, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 68984, 0, 3, 59954, 60374, 39944, 64934,
                                                 23144, 23414, 44264, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 69659, 0, 3, 61214, 61634, 41024, 66014,
                                                 23954, 24224, 45614, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 70334, 0, 3, 61634, 62054, 41384, 66554,
                                                 24224, 24494, 46064, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 71009, 0, 3, 62054, 62474, 41744, 67094,
                                                 24494, 24764, 46514, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71684, 3, 25304, 25314, 46979, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71705, 3, 25314, 25324, 46994, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71726, 3, 25324, 25334, 47009, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71747, 3, 25334, 25344, 47024, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71768, 3, 25344, 25354, 47039, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71789, 3, 25354, 25364, 47054, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71810, 3, 25364, 25374, 47069, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71831, 3, 25374, 25384, 47084, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71852, 3, 25384, 25394, 47099, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71873, 3, 25414, 25424, 47129, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71894, 3, 25424, 25434, 47144, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71915, 3, 25434, 25444, 47159, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71936, 3, 25444, 25454, 47174, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71957, 3, 25454, 25464, 47189, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71978, 3, 25464, 25474, 47204, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 71999, 3, 25474, 25484, 47219, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 72020, 3, 25484, 25494, 47234, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 72041, 3, 25494, 25504, 47249, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 72062, 0, 3, 46964, 71684, 47309, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72125, 0, 3, 46979, 71705, 47354, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72188, 0, 3, 46994, 71726, 47399, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72251, 0, 3, 47009, 71747, 47444, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72314, 0, 3, 47024, 71768, 47489, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72377, 0, 3, 47039, 71789, 47534, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72440, 0, 3, 47054, 71810, 47579, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72503, 0, 3, 47069, 71831, 47624, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72566, 0, 3, 47084, 71852, 47669, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72629, 0, 3, 47114, 71873, 47759, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72692, 0, 3, 47129, 71894, 47804, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72755, 0, 3, 47144, 71915, 47849, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72818, 0, 3, 47159, 71936, 47894, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72881, 0, 3, 47174, 71957, 47939, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 72944, 0, 3, 47189, 71978, 47984, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 73007, 0, 3, 47204, 71999, 48029, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 73070, 0, 3, 47219, 72020, 48074, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 73133, 0, 3, 47234, 72041, 48119, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 73196, 0, 3, 47264, 72062, 26184, 26244,
                                                 48164, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73322, 0, 3, 47309, 72125, 26244, 26304,
                                                 48254, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73448, 0, 3, 47354, 72188, 26304, 26364,
                                                 48344, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73574, 0, 3, 47399, 72251, 26364, 26424,
                                                 48434, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73700, 0, 3, 47444, 72314, 26424, 26484,
                                                 48524, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73826, 0, 3, 47489, 72377, 26484, 26544,
                                                 48614, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 73952, 0, 3, 47534, 72440, 26544, 26604,
                                                 48704, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74078, 0, 3, 47579, 72503, 26604, 26664,
                                                 48794, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74204, 0, 3, 47624, 72566, 26664, 26724,
                                                 48884, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74330, 0, 3, 47714, 72629, 26844, 26904,
                                                 48974, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74456, 0, 3, 47759, 72692, 26904, 26964,
                                                 49064, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74582, 0, 3, 47804, 72755, 26964, 27024,
                                                 49154, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74708, 0, 3, 47849, 72818, 27024, 27084,
                                                 49244, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74834, 0, 3, 47894, 72881, 27084, 27144,
                                                 49334, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 74960, 0, 3, 47939, 72944, 27144, 27204,
                                                 49424, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 75086, 0, 3, 47984, 73007, 27204, 27264,
                                                 49514, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 75212, 0, 3, 48029, 73070, 27264, 27324,
                                                 49604, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 75338, 0, 3, 48074, 73133, 27324, 27384,
                                                 49694, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 75464, 0, 3, 48254, 73448, 27504, 27604,
                                                 49934, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 75674, 0, 3, 48344, 73574, 27604, 27704,
                                                 50084, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 75884, 0, 3, 48434, 73700, 27704, 27804,
                                                 50234, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 76094, 0, 3, 48524, 73826, 27804, 27904,
                                                 50384, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 76304, 0, 3, 48614, 73952, 27904, 28004,
                                                 50534, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 76514, 0, 3, 48704, 74078, 28004, 28104,
                                                 50684, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 76724, 0, 3, 48794, 74204, 28104, 28204,
                                                 50834, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 76934, 0, 3, 49064, 74582, 28404, 28504,
                                                 51134, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 77144, 0, 3, 49154, 74708, 28504, 28604,
                                                 51284, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 77354, 0, 3, 49244, 74834, 28604, 28704,
                                                 51434, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 77564, 0, 3, 49334, 74960, 28704, 28804,
                                                 51584, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 77774, 0, 3, 49424, 75086, 28804, 28904,
                                                 51734, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 77984, 0, 3, 49514, 75212, 28904, 29004,
                                                 51884, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 78194, 0, 3, 49604, 75338, 29004, 29104,
                                                 52034, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 78404, 0, 3, 73196, 73322, 49784, 75464,
                                                 29304, 29454, 52184, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 78719, 0, 3, 73322, 73448, 49934, 75674,
                                                 29454, 29604, 52409, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 79034, 0, 3, 73448, 73574, 50084, 75884,
                                                 29604, 29754, 52634, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 79349, 0, 3, 73574, 73700, 50234, 76094,
                                                 29754, 29904, 52859, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 79664, 0, 3, 73700, 73826, 50384, 76304,
                                                 29904, 30054, 53084, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 79979, 0, 3, 73826, 73952, 50534, 76514,
                                                 30054, 30204, 53309, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 80294, 0, 3, 73952, 74078, 50684, 76724,
                                                 30204, 30354, 53534, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 80609, 0, 3, 74330, 74456, 50984, 76934,
                                                 30654, 30804, 53759, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 80924, 0, 3, 74456, 74582, 51134, 77144,
                                                 30804, 30954, 53984, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 81239, 0, 3, 74582, 74708, 51284, 77354,
                                                 30954, 31104, 54209, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 81554, 0, 3, 74708, 74834, 51434, 77564,
                                                 31104, 31254, 54434, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 81869, 0, 3, 74834, 74960, 51584, 77774,
                                                 31254, 31404, 54659, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 82184, 0, 3, 74960, 75086, 51734, 77984,
                                                 31404, 31554, 54884, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 82499, 0, 3, 75086, 75212, 51884, 78194,
                                                 31554, 31704, 55109, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 82814, 0, 3, 75464, 75674, 52409, 79034,
                                                 32004, 32214, 55649, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 83255, 0, 3, 75674, 75884, 52634, 79349,
                                                 32214, 32424, 55964, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 83696, 0, 3, 75884, 76094, 52859, 79664,
                                                 32424, 32634, 56279, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 84137, 0, 3, 76094, 76304, 53084, 79979,
                                                 32634, 32844, 56594, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 84578, 0, 3, 76304, 76514, 53309, 80294,
                                                 32844, 33054, 56909, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 85019, 0, 3, 76934, 77144, 53984, 81239,
                                                 33474, 33684, 57539, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 85460, 0, 3, 77144, 77354, 54209, 81554,
                                                 33684, 33894, 57854, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 85901, 0, 3, 77354, 77564, 54434, 81869,
                                                 33894, 34104, 58169, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 86342, 0, 3, 77564, 77774, 54659, 82184,
                                                 34104, 34314, 58484, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 86783, 0, 3, 77774, 77984, 54884, 82499,
                                                 34314, 34524, 58799, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 87224, 0, 3, 78404, 78719, 55334, 82814,
                                                 34944, 35224, 59114, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 87812, 0, 3, 78719, 79034, 55649, 83255,
                                                 35224, 35504, 59534, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 88400, 0, 3, 79034, 79349, 55964, 83696,
                                                 35504, 35784, 59954, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 88988, 0, 3, 79349, 79664, 56279, 84137,
                                                 35784, 36064, 60374, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 89576, 0, 3, 79664, 79979, 56594, 84578,
                                                 36064, 36344, 60794, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 90164, 0, 3, 80609, 80924, 57224, 85019,
                                                 36904, 37184, 61214, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 90752, 0, 3, 80924, 81239, 57539, 85460,
                                                 37184, 37464, 61634, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 91340, 0, 3, 81239, 81554, 57854, 85901,
                                                 37464, 37744, 62054, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 91928, 0, 3, 81554, 81869, 58169, 86342,
                                                 37744, 38024, 62474, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 92516, 0, 3, 81869, 82184, 58484, 86783,
                                                 38024, 38304, 62894, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 93104, 0, 3, 82814, 83255, 59534, 88400,
                                                 38864, 39224, 63854, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 93860, 0, 3, 83255, 83696, 59954, 88988,
                                                 39224, 39584, 64394, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 94616, 0, 3, 83696, 84137, 60374, 89576,
                                                 39584, 39944, 64934, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 95372, 0, 3, 85019, 85460, 61634, 91340,
                                                 40664, 41024, 66014, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 96128, 0, 3, 85460, 85901, 62054, 91928,
                                                 41024, 41384, 66554, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 96884, 0, 3, 85901, 86342, 62474, 92516,
                                                 41384, 41744, 67094, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 97640, 0, 3, 87224, 87812, 63314, 93104,
                                                 42464, 42914, 67634, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 98585, 0, 3, 87812, 88400, 63854, 93860,
                                                 42914, 43364, 68309, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 99530, 0, 3, 88400, 88988, 64394, 94616,
                                                 43364, 43814, 68984, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 100475, 0, 3, 90164, 90752, 65474,
                                                 95372, 44714, 45164, 69659, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 101420, 0, 3, 90752, 91340, 66014,
                                                 96128, 45164, 45614, 70334, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 102365, 0, 3, 91340, 91928, 66554,
                                                 96884, 45614, 46064, 71009, ncols, alpha, beta,
                                                 p);

            compute_prim_si_electron_repulsion_0(buffer, 103310, 3, 46964, 46979, 71705, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103338, 3, 46979, 46994, 71726, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103366, 3, 46994, 47009, 71747, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103394, 3, 47009, 47024, 71768, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103422, 3, 47024, 47039, 71789, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103450, 3, 47039, 47054, 71810, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103478, 3, 47054, 47069, 71831, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103506, 3, 47069, 47084, 71852, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103534, 3, 47114, 47129, 71894, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103562, 3, 47129, 47144, 71915, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103590, 3, 47144, 47159, 71936, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103618, 3, 47159, 47174, 71957, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103646, 3, 47174, 47189, 71978, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103674, 3, 47189, 47204, 71999, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103702, 3, 47204, 47219, 72020, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 103730, 3, 47219, 47234, 72041, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 103758, 0, 3, 71684, 103310, 72125,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 103842, 0, 3, 71705, 103338, 72188,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 103926, 0, 3, 71726, 103366, 72251,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104010, 0, 3, 71747, 103394, 72314,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104094, 0, 3, 71768, 103422, 72377,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104178, 0, 3, 71789, 103450, 72440,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104262, 0, 3, 71810, 103478, 72503,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104346, 0, 3, 71831, 103506, 72566,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104430, 0, 3, 71873, 103534, 72692,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104514, 0, 3, 71894, 103562, 72755,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104598, 0, 3, 71915, 103590, 72818,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104682, 0, 3, 71936, 103618, 72881,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104766, 0, 3, 71957, 103646, 72944,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104850, 0, 3, 71978, 103674, 73007,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 104934, 0, 3, 71999, 103702, 73070,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 105018, 0, 3, 72020, 103730, 73133,
                                                 ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 105102, 0, 3, 72125, 103842, 48164,
                                                 48254, 73448, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 105270, 0, 3, 72188, 103926, 48254,
                                                 48344, 73574, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 105438, 0, 3, 72251, 104010, 48344,
                                                 48434, 73700, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 105606, 0, 3, 72314, 104094, 48434,
                                                 48524, 73826, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 105774, 0, 3, 72377, 104178, 48524,
                                                 48614, 73952, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 105942, 0, 3, 72440, 104262, 48614,
                                                 48704, 74078, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106110, 0, 3, 72503, 104346, 48704,
                                                 48794, 74204, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106278, 0, 3, 72692, 104514, 48974,
                                                 49064, 74582, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106446, 0, 3, 72755, 104598, 49064,
                                                 49154, 74708, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106614, 0, 3, 72818, 104682, 49154,
                                                 49244, 74834, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106782, 0, 3, 72881, 104766, 49244,
                                                 49334, 74960, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 106950, 0, 3, 72944, 104850, 49334,
                                                 49424, 75086, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 107118, 0, 3, 73007, 104934, 49424,
                                                 49514, 75212, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 107286, 0, 3, 73070, 105018, 49514,
                                                 49604, 75338, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 107454, 0, 3, 73448, 105270, 49784,
                                                 49934, 75674, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 107734, 0, 3, 73574, 105438, 49934,
                                                 50084, 75884, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 108014, 0, 3, 73700, 105606, 50084,
                                                 50234, 76094, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 108294, 0, 3, 73826, 105774, 50234,
                                                 50384, 76304, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 108574, 0, 3, 73952, 105942, 50384,
                                                 50534, 76514, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 108854, 0, 3, 74078, 106110, 50534,
                                                 50684, 76724, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 109134, 0, 3, 74582, 106446, 50984,
                                                 51134, 77144, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 109414, 0, 3, 74708, 106614, 51134,
                                                 51284, 77354, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 109694, 0, 3, 74834, 106782, 51284,
                                                 51434, 77564, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 109974, 0, 3, 74960, 106950, 51434,
                                                 51584, 77774, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 110254, 0, 3, 75086, 107118, 51584,
                                                 51734, 77984, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 110534, 0, 3, 75212, 107286, 51734,
                                                 51884, 78194, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 110814, 0, 3, 105102, 105270, 75674,
                                                 107734, 52184, 52409, 79034, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 111234, 0, 3, 105270, 105438, 75884,
                                                 108014, 52409, 52634, 79349, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 111654, 0, 3, 105438, 105606, 76094,
                                                 108294, 52634, 52859, 79664, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 112074, 0, 3, 105606, 105774, 76304,
                                                 108574, 52859, 53084, 79979, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 112494, 0, 3, 105774, 105942, 76514,
                                                 108854, 53084, 53309, 80294, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 112914, 0, 3, 106278, 106446, 77144,
                                                 109414, 53759, 53984, 81239, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 113334, 0, 3, 106446, 106614, 77354,
                                                 109694, 53984, 54209, 81554, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 113754, 0, 3, 106614, 106782, 77564,
                                                 109974, 54209, 54434, 81869, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 114174, 0, 3, 106782, 106950, 77774,
                                                 110254, 54434, 54659, 82184, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 114594, 0, 3, 106950, 107118, 77984,
                                                 110534, 54659, 54884, 82499, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 115014, 0, 3, 107454, 107734, 79034,
                                                 111234, 55334, 55649, 83255, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 115602, 0, 3, 107734, 108014, 79349,
                                                 111654, 55649, 55964, 83696, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 116190, 0, 3, 108014, 108294, 79664,
                                                 112074, 55964, 56279, 84137, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 116778, 0, 3, 108294, 108574, 79979,
                                                 112494, 56279, 56594, 84578, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 117366, 0, 3, 109134, 109414, 81239,
                                                 113334, 57224, 57539, 85460, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 117954, 0, 3, 109414, 109694, 81554,
                                                 113754, 57539, 57854, 85901, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 118542, 0, 3, 109694, 109974, 81869,
                                                 114174, 57854, 58169, 86342, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 119130, 0, 3, 109974, 110254, 82184,
                                                 114594, 58169, 58484, 86783, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 119718, 0, 3, 110814, 111234, 83255,
                                                 115602, 59114, 59534, 88400, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 120502, 0, 3, 111234, 111654, 83696,
                                                 116190, 59534, 59954, 88988, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 121286, 0, 3, 111654, 112074, 84137,
                                                 116778, 59954, 60374, 89576, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 122070, 0, 3, 112914, 113334, 85460,
                                                 117954, 61214, 61634, 91340, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 122854, 0, 3, 113334, 113754, 85901,
                                                 118542, 61634, 62054, 91928, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 123638, 0, 3, 113754, 114174, 86342,
                                                 119130, 62054, 62474, 92516, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 124422, 0, 3, 115014, 115602, 88400,
                                                 120502, 63314, 63854, 93860, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 125430, 0, 3, 115602, 116190, 88988,
                                                 121286, 63854, 64394, 94616, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 126438, 0, 3, 117366, 117954, 91340,
                                                 122854, 65474, 66014, 96128, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 127446, 0, 3, 117954, 118542, 91928,
                                                 123638, 66014, 66554, 96884, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 128454, 0, 3, 119718, 120502, 93860,
                                                 125430, 67634, 68309, 99530, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 129714, 0, 3, 122070, 122854, 96128,
                                                 127446, 69659, 70334, 102365, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 130974, 3, 71684, 71705, 103338, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131010, 3, 71705, 71726, 103366, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131046, 3, 71726, 71747, 103394, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131082, 3, 71747, 71768, 103422, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131118, 3, 71768, 71789, 103450, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131154, 3, 71789, 71810, 103478, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131190, 3, 71810, 71831, 103506, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131226, 3, 71873, 71894, 103562, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131262, 3, 71894, 71915, 103590, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131298, 3, 71915, 71936, 103618, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131334, 3, 71936, 71957, 103646, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131370, 3, 71957, 71978, 103674, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131406, 3, 71978, 71999, 103702, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 131442, 3, 71999, 72020, 103730, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 131478, 0, 3, 103310, 130974, 103842,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 131586, 0, 3, 103338, 131010, 103926,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 131694, 0, 3, 103366, 131046, 104010,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 131802, 0, 3, 103394, 131082, 104094,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 131910, 0, 3, 103422, 131118, 104178,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132018, 0, 3, 103450, 131154, 104262,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132126, 0, 3, 103478, 131190, 104346,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132234, 0, 3, 103534, 131226, 104514,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132342, 0, 3, 103562, 131262, 104598,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132450, 0, 3, 103590, 131298, 104682,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132558, 0, 3, 103618, 131334, 104766,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132666, 0, 3, 103646, 131370, 104850,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132774, 0, 3, 103674, 131406, 104934,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 132882, 0, 3, 103702, 131442, 105018,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 132990, 0, 3, 103758, 131478, 73196,
                                                 73322, 105102, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 133206, 0, 3, 103842, 131586, 73322,
                                                 73448, 105270, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 133422, 0, 3, 103926, 131694, 73448,
                                                 73574, 105438, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 133638, 0, 3, 104010, 131802, 73574,
                                                 73700, 105606, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 133854, 0, 3, 104094, 131910, 73700,
                                                 73826, 105774, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 134070, 0, 3, 104178, 132018, 73826,
                                                 73952, 105942, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 134286, 0, 3, 104262, 132126, 73952,
                                                 74078, 106110, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 134502, 0, 3, 104430, 132234, 74330,
                                                 74456, 106278, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 134718, 0, 3, 104514, 132342, 74456,
                                                 74582, 106446, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 134934, 0, 3, 104598, 132450, 74582,
                                                 74708, 106614, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 135150, 0, 3, 104682, 132558, 74708,
                                                 74834, 106782, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 135366, 0, 3, 104766, 132666, 74834,
                                                 74960, 106950, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 135582, 0, 3, 104850, 132774, 74960,
                                                 75086, 107118, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 135798, 0, 3, 104934, 132882, 75086,
                                                 75212, 107286, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 136014, 0, 3, 105270, 133422, 75464,
                                                 75674, 107734, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 136374, 0, 3, 105438, 133638, 75674,
                                                 75884, 108014, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 136734, 0, 3, 105606, 133854, 75884,
                                                 76094, 108294, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 137094, 0, 3, 105774, 134070, 76094,
                                                 76304, 108574, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 137454, 0, 3, 105942, 134286, 76304,
                                                 76514, 108854, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 137814, 0, 3, 106446, 134934, 76934,
                                                 77144, 109414, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 138174, 0, 3, 106614, 135150, 77144,
                                                 77354, 109694, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 138534, 0, 3, 106782, 135366, 77354,
                                                 77564, 109974, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 138894, 0, 3, 106950, 135582, 77564,
                                                 77774, 110254, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 139254, 0, 3, 107118, 135798, 77774,
                                                 77984, 110534, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 139614, 0, 3, 132990, 133206, 107454,
                                                 136014, 78404, 78719, 110814, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 140154, 0, 3, 133206, 133422, 107734,
                                                 136374, 78719, 79034, 111234, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 140694, 0, 3, 133422, 133638, 108014,
                                                 136734, 79034, 79349, 111654, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 141234, 0, 3, 133638, 133854, 108294,
                                                 137094, 79349, 79664, 112074, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 141774, 0, 3, 133854, 134070, 108574,
                                                 137454, 79664, 79979, 112494, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 142314, 0, 3, 134502, 134718, 109134,
                                                 137814, 80609, 80924, 112914, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 142854, 0, 3, 134718, 134934, 109414,
                                                 138174, 80924, 81239, 113334, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 143394, 0, 3, 134934, 135150, 109694,
                                                 138534, 81239, 81554, 113754, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 143934, 0, 3, 135150, 135366, 109974,
                                                 138894, 81554, 81869, 114174, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 144474, 0, 3, 135366, 135582, 110254,
                                                 139254, 81869, 82184, 114594, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 145014, 0, 3, 136014, 136374, 111234,
                                                 140694, 82814, 83255, 115602, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 145770, 0, 3, 136374, 136734, 111654,
                                                 141234, 83255, 83696, 116190, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 146526, 0, 3, 136734, 137094, 112074,
                                                 141774, 83696, 84137, 116778, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 147282, 0, 3, 137814, 138174, 113334,
                                                 143394, 85019, 85460, 117954, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 148038, 0, 3, 138174, 138534, 113754,
                                                 143934, 85460, 85901, 118542, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 148794, 0, 3, 138534, 138894, 114174,
                                                 144474, 85901, 86342, 119130, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 149550, 0, 3, 139614, 140154, 115014,
                                                 145014, 87224, 87812, 119718, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 150558, 0, 3, 140154, 140694, 115602,
                                                 145770, 87812, 88400, 120502, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 151566, 0, 3, 140694, 141234, 116190,
                                                 146526, 88400, 88988, 121286, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 152574, 0, 3, 142314, 142854, 117366,
                                                 147282, 90164, 90752, 122070, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 153582, 0, 3, 142854, 143394, 117954,
                                                 148038, 90752, 91340, 122854, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 154590, 0, 3, 143394, 143934, 118542,
                                                 148794, 91340, 91928, 123638, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 155598, 0, 3, 145014, 145770, 120502,
                                                 151566, 93104, 93860, 125430, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 156894, 0, 3, 147282, 148038, 122854,
                                                 154590, 95372, 96128, 127446, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 158190, 0, 3, 149550, 150558, 124422,
                                                 155598, 97640, 98585, 128454, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 159810, 0, 3, 152574, 153582, 126438,
                                                 156894, 100475, 101420, 129714, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_k_x(buffer, 161430, 152574, 159810, 1, 36, ncols, alpha);

            simdgeo::geom_k_y(buffer, 162726, 152574, 159810, 1, 36, ncols, alpha);

            simdgeo::geom_k_z(buffer, 164022, 152574, 159810, 1, 36, ncols, alpha);

            simdgeo::geom_k_x(buffer, 165318, 149550, 158190, 1, 36, ncols, alpha);

            simdgeo::geom_k_y(buffer, 166614, 149550, 158190, 1, 36, ncols, alpha);

            simdgeo::geom_k_z(buffer, 167910, 149550, 158190, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 169206, 165318, 3888, ncols);

            simdfunc::contract_primitives(buffer, 173094, 161430, 3888, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 176982, 173094, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 176982, 15, nmax);

    simdtrf::transform_k_inner(buffer, 176982, 174390, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 225 * nvalues, nvalues, buffer, 176982, 15, nmax);

    simdtrf::transform_k_inner(buffer, 176982, 175686, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 450 * nvalues, nvalues, buffer, 176982, 15, nmax);

    simdtrf::transform_k_inner(buffer, 176982, 169206, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 675 * nvalues, nvalues, buffer, 176982, 15, nmax);

    simdtrf::transform_k_inner(buffer, 176982, 170502, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 900 * nvalues, nvalues, buffer, 176982, 15, nmax);

    simdtrf::transform_k_inner(buffer, 176982, 171798, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 1125 * nvalues, nvalues, buffer, 176982, 15, nmax);
}

}  // namespace simdt2ceri
