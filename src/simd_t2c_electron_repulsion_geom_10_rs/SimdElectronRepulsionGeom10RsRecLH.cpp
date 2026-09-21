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


#include "SimdElectronRepulsionGeom10RsRecLH.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
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
#include "SimdGeometryL1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_lh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_lh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 105829, 99664, 5670, nvalues);

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
                                                9, 10, 11, 12, 13, 14}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 21, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 96, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 99, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 102, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 105, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 108, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 111, 0, 35, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 7, 8, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 8, 9, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 9, 10, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 10, 11, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 11, 12, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 12, 13, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 13, 14, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 14, 15, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 15, 16, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 16, 17, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 17, 18, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 18, 19, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 22, 23, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 23, 24, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 24, 25, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 204, 0, 25, 26, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 210, 0, 26, 27, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 216, 0, 27, 28, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 222, 0, 28, 29, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 228, 0, 29, 30, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 234, 0, 30, 31, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 240, 0, 31, 32, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 246, 0, 32, 33, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 252, 0, 33, 34, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 36, 39, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 39, 42, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 42, 45, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 45, 48, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 48, 51, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 51, 54, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 54, 57, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 57, 60, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 60, 63, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 63, 66, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 66, 69, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 75, 78, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 78, 81, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 81, 84, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 398, 0, 84, 87, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 408, 0, 87, 90, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 418, 0, 90, 93, 222, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 428, 0, 93, 96, 228, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 438, 0, 96, 99, 234, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 448, 0, 99, 102, 240, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 458, 0, 102, 105, 246, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 468, 0, 105, 108, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 478, 0, 114, 120, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 493, 0, 120, 126, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 508, 0, 126, 132, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 523, 0, 132, 138, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 538, 0, 138, 144, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 553, 0, 144, 150, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 568, 0, 150, 156, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 583, 0, 156, 162, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 598, 0, 162, 168, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 613, 0, 168, 174, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 628, 0, 186, 192, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 643, 0, 192, 198, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 658, 0, 198, 204, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 673, 0, 204, 210, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 688, 0, 210, 216, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 703, 0, 216, 222, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 718, 0, 222, 228, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 733, 0, 228, 234, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 748, 0, 234, 240, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 763, 0, 240, 246, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 778, 0, 258, 268, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 799, 0, 268, 278, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 820, 0, 278, 288, 523, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 841, 0, 288, 298, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 862, 0, 298, 308, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 883, 0, 308, 318, 568, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 904, 0, 318, 328, 583, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 925, 0, 328, 338, 598, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 946, 0, 338, 348, 613, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 967, 0, 368, 378, 643, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 988, 0, 378, 388, 658, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1009, 0, 388, 398, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1030, 0, 398, 408, 688, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1051, 0, 408, 418, 703, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1072, 0, 418, 428, 718, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1093, 0, 428, 438, 733, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1114, 0, 438, 448, 748, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1135, 0, 448, 458, 763, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1156, 0, 478, 493, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1184, 0, 493, 508, 820, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1212, 0, 508, 523, 841, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1240, 0, 523, 538, 862, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1268, 0, 538, 553, 883, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1296, 0, 553, 568, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1324, 0, 568, 583, 925, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1352, 0, 583, 598, 946, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1380, 0, 628, 643, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1408, 0, 643, 658, 1009, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1436, 0, 658, 673, 1030, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1464, 0, 673, 688, 1051, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1492, 0, 688, 703, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1520, 0, 703, 718, 1093, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1548, 0, 718, 733, 1114, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1576, 0, 733, 748, 1135, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1604, 0, 778, 799, 1184, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1640, 0, 799, 820, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1676, 0, 820, 841, 1240, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1712, 0, 841, 862, 1268, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1748, 0, 862, 883, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1784, 0, 883, 904, 1324, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1820, 0, 904, 925, 1352, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1856, 0, 967, 988, 1408, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1892, 0, 988, 1009, 1436, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1928, 0, 1009, 1030, 1464, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1964, 0, 1030, 1051, 1492, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2000, 0, 1051, 1072, 1520, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2036, 0, 1072, 1093, 1548, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2072, 0, 1093, 1114, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2108, 0, 1156, 1184, 1640, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2153, 0, 1184, 1212, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2198, 0, 1212, 1240, 1712, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2243, 0, 1240, 1268, 1748, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2288, 0, 1268, 1296, 1784, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2333, 0, 1296, 1324, 1820, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2378, 0, 1380, 1408, 1892, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2423, 0, 1408, 1436, 1928, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2468, 0, 1436, 1464, 1964, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2513, 0, 1464, 1492, 2000, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2558, 0, 1492, 1520, 2036, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2603, 0, 1520, 1548, 2072, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2648, 0, 1604, 1640, 2153, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2703, 0, 1640, 1676, 2198, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2758, 0, 1676, 1712, 2243, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2813, 0, 1712, 1748, 2288, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2868, 0, 1748, 1784, 2333, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2923, 0, 1856, 1892, 2423, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2978, 0, 1892, 1928, 2468, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3033, 0, 1928, 1964, 2513, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3088, 0, 1964, 2000, 2558, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3143, 0, 2000, 2036, 2603, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 3198, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3201, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3204, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3207, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3210, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3213, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3216, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3219, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3222, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3225, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3228, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3231, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3234, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3237, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3240, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3243, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3246, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3249, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3252, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3255, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3258, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3261, 3, 35, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 3264, 3, 9, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3273, 3, 10, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3282, 3, 11, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3291, 3, 12, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3300, 3, 13, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3309, 3, 14, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3318, 3, 15, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3327, 3, 16, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3336, 3, 17, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3345, 3, 18, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3354, 3, 19, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3363, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3372, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3381, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3390, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3399, 3, 28, 93, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3408, 3, 29, 96, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3417, 3, 30, 99, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3426, 3, 31, 102, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3435, 3, 32, 105, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3444, 3, 33, 108, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3453, 3, 34, 111, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3462, 0, 3, 39, 3264, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3480, 0, 3, 42, 3273, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3498, 0, 3, 45, 3282, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3516, 0, 3, 48, 3291, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3534, 0, 3, 51, 3300, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3552, 0, 3, 54, 3309, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3570, 0, 3, 57, 3318, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3588, 0, 3, 60, 3327, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3606, 0, 3, 63, 3336, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3624, 0, 3, 66, 3345, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3642, 0, 3, 69, 3354, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3660, 0, 3, 78, 3363, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3678, 0, 3, 81, 3372, 198, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3696, 0, 3, 84, 3381, 204, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3714, 0, 3, 87, 3390, 210, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3732, 0, 3, 90, 3399, 216, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3750, 0, 3, 93, 3408, 222, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3768, 0, 3, 96, 3417, 228, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3786, 0, 3, 99, 3426, 234, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3804, 0, 3, 102, 3435, 240, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3822, 0, 3, 105, 3444, 246, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3840, 0, 3, 108, 3453, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3858, 0, 3, 114, 3462, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3888, 0, 3, 120, 3480, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3918, 0, 3, 126, 3498, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3948, 0, 3, 132, 3516, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3978, 0, 3, 138, 3534, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4008, 0, 3, 144, 3552, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4038, 0, 3, 150, 3570, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4068, 0, 3, 156, 3588, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4098, 0, 3, 162, 3606, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4128, 0, 3, 168, 3624, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4158, 0, 3, 174, 3642, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4188, 0, 3, 186, 3660, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4218, 0, 3, 192, 3678, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4248, 0, 3, 198, 3696, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4278, 0, 3, 204, 3714, 398, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4308, 0, 3, 210, 3732, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4338, 0, 3, 216, 3750, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4368, 0, 3, 222, 3768, 428, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4398, 0, 3, 228, 3786, 438, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4428, 0, 3, 234, 3804, 448, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4458, 0, 3, 240, 3822, 458, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4488, 0, 3, 246, 3840, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4518, 0, 3, 268, 3918, 493, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4563, 0, 3, 278, 3948, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4608, 0, 3, 288, 3978, 523, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4653, 0, 3, 298, 4008, 538, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4698, 0, 3, 308, 4038, 553, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4743, 0, 3, 318, 4068, 568, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4788, 0, 3, 328, 4098, 583, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4833, 0, 3, 338, 4128, 598, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4878, 0, 3, 348, 4158, 613, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4923, 0, 3, 378, 4248, 643, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4968, 0, 3, 388, 4278, 658, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5013, 0, 3, 398, 4308, 673, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5058, 0, 3, 408, 4338, 688, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5103, 0, 3, 418, 4368, 703, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5148, 0, 3, 428, 4398, 718, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5193, 0, 3, 438, 4428, 733, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5238, 0, 3, 448, 4458, 748, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5283, 0, 3, 458, 4488, 763, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5328, 0, 3, 478, 4518, 778, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5391, 0, 3, 493, 4563, 799, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5454, 0, 3, 508, 4608, 820, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5517, 0, 3, 523, 4653, 841, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5580, 0, 3, 538, 4698, 862, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5643, 0, 3, 553, 4743, 883, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5706, 0, 3, 568, 4788, 904, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5769, 0, 3, 583, 4833, 925, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5832, 0, 3, 598, 4878, 946, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5895, 0, 3, 628, 4923, 967, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5958, 0, 3, 643, 4968, 988, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6021, 0, 3, 658, 5013, 1009, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6084, 0, 3, 673, 5058, 1030, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6147, 0, 3, 688, 5103, 1051, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6210, 0, 3, 703, 5148, 1072, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6273, 0, 3, 718, 5193, 1093, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6336, 0, 3, 733, 5238, 1114, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6399, 0, 3, 748, 5283, 1135, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6462, 0, 3, 799, 5454, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6546, 0, 3, 820, 5517, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6630, 0, 3, 841, 5580, 1240, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6714, 0, 3, 862, 5643, 1268, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6798, 0, 3, 883, 5706, 1296, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6882, 0, 3, 904, 5769, 1324, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6966, 0, 3, 925, 5832, 1352, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7050, 0, 3, 988, 6021, 1408, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7134, 0, 3, 1009, 6084, 1436, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7218, 0, 3, 1030, 6147, 1464, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7302, 0, 3, 1051, 6210, 1492, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7386, 0, 3, 1072, 6273, 1520, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7470, 0, 3, 1093, 6336, 1548, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7554, 0, 3, 1114, 6399, 1576, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7638, 0, 3, 1156, 6462, 1604, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7746, 0, 3, 1184, 6546, 1640, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7854, 0, 3, 1212, 6630, 1676, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7962, 0, 3, 1240, 6714, 1712, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8070, 0, 3, 1268, 6798, 1748, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8178, 0, 3, 1296, 6882, 1784, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8286, 0, 3, 1324, 6966, 1820, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8394, 0, 3, 1380, 7050, 1856, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8502, 0, 3, 1408, 7134, 1892, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8610, 0, 3, 1436, 7218, 1928, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8718, 0, 3, 1464, 7302, 1964, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8826, 0, 3, 1492, 7386, 2000, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8934, 0, 3, 1520, 7470, 2036, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9042, 0, 3, 1548, 7554, 2072, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9150, 0, 3, 1640, 7854, 2153, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9285, 0, 3, 1676, 7962, 2198, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9420, 0, 3, 1712, 8070, 2243, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9555, 0, 3, 1748, 8178, 2288, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9690, 0, 3, 1784, 8286, 2333, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9825, 0, 3, 1892, 8610, 2423, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9960, 0, 3, 1928, 8718, 2468, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10095, 0, 3, 1964, 8826, 2513, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10230, 0, 3, 2000, 8934, 2558, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10365, 0, 3, 2036, 9042, 2603, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 10500, 0, 3, 2108, 9150, 2648, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 10665, 0, 3, 2153, 9285, 2703, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 10830, 0, 3, 2198, 9420, 2758, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 10995, 0, 3, 2243, 9555, 2813, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11160, 0, 3, 2288, 9690, 2868, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11325, 0, 3, 2378, 9825, 2923, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11490, 0, 3, 2423, 9960, 2978, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11655, 0, 3, 2468, 10095, 3033, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11820, 0, 3, 2513, 10230, 3088, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11985, 0, 3, 2558, 10365, 3143, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 12150, 3, 9, 10, 3201, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12156, 3, 10, 11, 3204, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12162, 3, 11, 12, 3207, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12168, 3, 12, 13, 3210, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12174, 3, 13, 14, 3213, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12180, 3, 14, 15, 3216, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12186, 3, 15, 16, 3219, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12192, 3, 16, 17, 3222, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12198, 3, 17, 18, 3225, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12204, 3, 18, 19, 3228, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12210, 3, 24, 25, 3234, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12216, 3, 25, 26, 3237, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12222, 3, 26, 27, 3240, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12228, 3, 27, 28, 3243, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12234, 3, 28, 29, 3246, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12240, 3, 29, 30, 3249, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12246, 3, 30, 31, 3252, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12252, 3, 31, 32, 3255, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12258, 3, 32, 33, 3258, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12264, 3, 33, 34, 3261, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 12270, 0, 3, 3198, 12150, 3273, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12288, 0, 3, 3201, 12156, 3282, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12306, 0, 3, 3204, 12162, 3291, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12324, 0, 3, 3207, 12168, 3300, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12342, 0, 3, 3210, 12174, 3309, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12360, 0, 3, 3213, 12180, 3318, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12378, 0, 3, 3216, 12186, 3327, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12396, 0, 3, 3219, 12192, 3336, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12414, 0, 3, 3222, 12198, 3345, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12432, 0, 3, 3225, 12204, 3354, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12450, 0, 3, 3231, 12210, 3372, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12468, 0, 3, 3234, 12216, 3381, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12486, 0, 3, 3237, 12222, 3390, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12504, 0, 3, 3240, 12228, 3399, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12522, 0, 3, 3243, 12234, 3408, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12540, 0, 3, 3246, 12240, 3417, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12558, 0, 3, 3249, 12246, 3426, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12576, 0, 3, 3252, 12252, 3435, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12594, 0, 3, 3255, 12258, 3444, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12612, 0, 3, 3258, 12264, 3453, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 12630, 0, 3, 3264, 12270, 114, 120,
                                                 3480, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12666, 0, 3, 3273, 12288, 120, 126,
                                                 3498, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12702, 0, 3, 3282, 12306, 126, 132,
                                                 3516, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12738, 0, 3, 3291, 12324, 132, 138,
                                                 3534, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12774, 0, 3, 3300, 12342, 138, 144,
                                                 3552, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12810, 0, 3, 3309, 12360, 144, 150,
                                                 3570, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12846, 0, 3, 3318, 12378, 150, 156,
                                                 3588, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12882, 0, 3, 3327, 12396, 156, 162,
                                                 3606, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12918, 0, 3, 3336, 12414, 162, 168,
                                                 3624, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12954, 0, 3, 3345, 12432, 168, 174,
                                                 3642, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 12990, 0, 3, 3363, 12450, 186, 192,
                                                 3678, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13026, 0, 3, 3372, 12468, 192, 198,
                                                 3696, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13062, 0, 3, 3381, 12486, 198, 204,
                                                 3714, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13098, 0, 3, 3390, 12504, 204, 210,
                                                 3732, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13134, 0, 3, 3399, 12522, 210, 216,
                                                 3750, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13170, 0, 3, 3408, 12540, 216, 222,
                                                 3768, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13206, 0, 3, 3417, 12558, 222, 228,
                                                 3786, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13242, 0, 3, 3426, 12576, 228, 234,
                                                 3804, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13278, 0, 3, 3435, 12594, 234, 240,
                                                 3822, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13314, 0, 3, 3444, 12612, 240, 246,
                                                 3840, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13350, 0, 3, 3480, 12666, 258, 268,
                                                 3918, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13410, 0, 3, 3498, 12702, 268, 278,
                                                 3948, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13470, 0, 3, 3516, 12738, 278, 288,
                                                 3978, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13530, 0, 3, 3534, 12774, 288, 298,
                                                 4008, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13590, 0, 3, 3552, 12810, 298, 308,
                                                 4038, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13650, 0, 3, 3570, 12846, 308, 318,
                                                 4068, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13710, 0, 3, 3588, 12882, 318, 328,
                                                 4098, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13770, 0, 3, 3606, 12918, 328, 338,
                                                 4128, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13830, 0, 3, 3624, 12954, 338, 348,
                                                 4158, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13890, 0, 3, 3678, 13026, 368, 378,
                                                 4248, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13950, 0, 3, 3696, 13062, 378, 388,
                                                 4278, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14010, 0, 3, 3714, 13098, 388, 398,
                                                 4308, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14070, 0, 3, 3732, 13134, 398, 408,
                                                 4338, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14130, 0, 3, 3750, 13170, 408, 418,
                                                 4368, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14190, 0, 3, 3768, 13206, 418, 428,
                                                 4398, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14250, 0, 3, 3786, 13242, 428, 438,
                                                 4428, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14310, 0, 3, 3804, 13278, 438, 448,
                                                 4458, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14370, 0, 3, 3822, 13314, 448, 458,
                                                 4488, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14430, 0, 3, 12630, 12666, 3918, 13410,
                                                 478, 493, 4563, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14520, 0, 3, 12666, 12702, 3948, 13470,
                                                 493, 508, 4608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14610, 0, 3, 12702, 12738, 3978, 13530,
                                                 508, 523, 4653, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14700, 0, 3, 12738, 12774, 4008, 13590,
                                                 523, 538, 4698, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14790, 0, 3, 12774, 12810, 4038, 13650,
                                                 538, 553, 4743, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14880, 0, 3, 12810, 12846, 4068, 13710,
                                                 553, 568, 4788, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 14970, 0, 3, 12846, 12882, 4098, 13770,
                                                 568, 583, 4833, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15060, 0, 3, 12882, 12918, 4128, 13830,
                                                 583, 598, 4878, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15150, 0, 3, 12990, 13026, 4248, 13950,
                                                 628, 643, 4968, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15240, 0, 3, 13026, 13062, 4278, 14010,
                                                 643, 658, 5013, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15330, 0, 3, 13062, 13098, 4308, 14070,
                                                 658, 673, 5058, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15420, 0, 3, 13098, 13134, 4338, 14130,
                                                 673, 688, 5103, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15510, 0, 3, 13134, 13170, 4368, 14190,
                                                 688, 703, 5148, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15600, 0, 3, 13170, 13206, 4398, 14250,
                                                 703, 718, 5193, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15690, 0, 3, 13206, 13242, 4428, 14310,
                                                 718, 733, 5238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15780, 0, 3, 13242, 13278, 4458, 14370,
                                                 733, 748, 5283, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15870, 0, 3, 13350, 13410, 4563, 14520,
                                                 778, 799, 5454, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15996, 0, 3, 13410, 13470, 4608, 14610,
                                                 799, 820, 5517, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16122, 0, 3, 13470, 13530, 4653, 14700,
                                                 820, 841, 5580, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16248, 0, 3, 13530, 13590, 4698, 14790,
                                                 841, 862, 5643, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16374, 0, 3, 13590, 13650, 4743, 14880,
                                                 862, 883, 5706, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16500, 0, 3, 13650, 13710, 4788, 14970,
                                                 883, 904, 5769, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16626, 0, 3, 13710, 13770, 4833, 15060,
                                                 904, 925, 5832, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16752, 0, 3, 13890, 13950, 4968, 15240,
                                                 967, 988, 6021, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 16878, 0, 3, 13950, 14010, 5013, 15330,
                                                 988, 1009, 6084, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17004, 0, 3, 14010, 14070, 5058, 15420,
                                                 1009, 1030, 6147, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17130, 0, 3, 14070, 14130, 5103, 15510,
                                                 1030, 1051, 6210, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17256, 0, 3, 14130, 14190, 5148, 15600,
                                                 1051, 1072, 6273, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17382, 0, 3, 14190, 14250, 5193, 15690,
                                                 1072, 1093, 6336, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17508, 0, 3, 14250, 14310, 5238, 15780,
                                                 1093, 1114, 6399, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17634, 0, 3, 14430, 14520, 5454, 15996,
                                                 1156, 1184, 6546, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17802, 0, 3, 14520, 14610, 5517, 16122,
                                                 1184, 1212, 6630, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17970, 0, 3, 14610, 14700, 5580, 16248,
                                                 1212, 1240, 6714, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18138, 0, 3, 14700, 14790, 5643, 16374,
                                                 1240, 1268, 6798, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18306, 0, 3, 14790, 14880, 5706, 16500,
                                                 1268, 1296, 6882, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18474, 0, 3, 14880, 14970, 5769, 16626,
                                                 1296, 1324, 6966, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18642, 0, 3, 15150, 15240, 6021, 16878,
                                                 1380, 1408, 7134, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18810, 0, 3, 15240, 15330, 6084, 17004,
                                                 1408, 1436, 7218, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 18978, 0, 3, 15330, 15420, 6147, 17130,
                                                 1436, 1464, 7302, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19146, 0, 3, 15420, 15510, 6210, 17256,
                                                 1464, 1492, 7386, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19314, 0, 3, 15510, 15600, 6273, 17382,
                                                 1492, 1520, 7470, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19482, 0, 3, 15600, 15690, 6336, 17508,
                                                 1520, 1548, 7554, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19650, 0, 3, 15870, 15996, 6546, 17802,
                                                 1604, 1640, 7854, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19866, 0, 3, 15996, 16122, 6630, 17970,
                                                 1640, 1676, 7962, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20082, 0, 3, 16122, 16248, 6714, 18138,
                                                 1676, 1712, 8070, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20298, 0, 3, 16248, 16374, 6798, 18306,
                                                 1712, 1748, 8178, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20514, 0, 3, 16374, 16500, 6882, 18474,
                                                 1748, 1784, 8286, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20730, 0, 3, 16752, 16878, 7134, 18810,
                                                 1856, 1892, 8610, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20946, 0, 3, 16878, 17004, 7218, 18978,
                                                 1892, 1928, 8718, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21162, 0, 3, 17004, 17130, 7302, 19146,
                                                 1928, 1964, 8826, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21378, 0, 3, 17130, 17256, 7386, 19314,
                                                 1964, 2000, 8934, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 21594, 0, 3, 17256, 17382, 7470, 19482,
                                                 2000, 2036, 9042, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 21810, 0, 3, 17634, 17802, 7854, 19866,
                                                 2108, 2153, 9285, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22080, 0, 3, 17802, 17970, 7962, 20082,
                                                 2153, 2198, 9420, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22350, 0, 3, 17970, 18138, 8070, 20298,
                                                 2198, 2243, 9555, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22620, 0, 3, 18138, 18306, 8178, 20514,
                                                 2243, 2288, 9690, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22890, 0, 3, 18642, 18810, 8610, 20946,
                                                 2378, 2423, 9960, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23160, 0, 3, 18810, 18978, 8718, 21162,
                                                 2423, 2468, 10095, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23430, 0, 3, 18978, 19146, 8826, 21378,
                                                 2468, 2513, 10230, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 23700, 0, 3, 19146, 19314, 8934, 21594,
                                                 2513, 2558, 10365, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 23970, 0, 3, 19650, 19866, 9285, 22080,
                                                 2648, 2703, 10830, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 24300, 0, 3, 19866, 20082, 9420, 22350,
                                                 2703, 2758, 10995, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 24630, 0, 3, 20082, 20298, 9555, 22620,
                                                 2758, 2813, 11160, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 24960, 0, 3, 20730, 20946, 9960, 23160,
                                                 2923, 2978, 11655, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 25290, 0, 3, 20946, 21162, 10095, 23430,
                                                 2978, 3033, 11820, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 25620, 0, 3, 21162, 21378, 10230, 23700,
                                                 3033, 3088, 11985, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25950, 3, 3198, 3201, 12156, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25960, 3, 3201, 3204, 12162, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25970, 3, 3204, 3207, 12168, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25980, 3, 3207, 3210, 12174, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 25990, 3, 3210, 3213, 12180, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26000, 3, 3213, 3216, 12186, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26010, 3, 3216, 3219, 12192, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26020, 3, 3219, 3222, 12198, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26030, 3, 3222, 3225, 12204, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26040, 3, 3231, 3234, 12216, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26050, 3, 3234, 3237, 12222, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26060, 3, 3237, 3240, 12228, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26070, 3, 3240, 3243, 12234, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26080, 3, 3243, 3246, 12240, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26090, 3, 3246, 3249, 12246, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26100, 3, 3249, 3252, 12252, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26110, 3, 3252, 3255, 12258, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 26120, 3, 3255, 3258, 12264, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 26130, 0, 3, 12150, 25950, 12288, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26160, 0, 3, 12156, 25960, 12306, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26190, 0, 3, 12162, 25970, 12324, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26220, 0, 3, 12168, 25980, 12342, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26250, 0, 3, 12174, 25990, 12360, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26280, 0, 3, 12180, 26000, 12378, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26310, 0, 3, 12186, 26010, 12396, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26340, 0, 3, 12192, 26020, 12414, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26370, 0, 3, 12198, 26030, 12432, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26400, 0, 3, 12210, 26040, 12468, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26430, 0, 3, 12216, 26050, 12486, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26460, 0, 3, 12222, 26060, 12504, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26490, 0, 3, 12228, 26070, 12522, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26520, 0, 3, 12234, 26080, 12540, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26550, 0, 3, 12240, 26090, 12558, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26580, 0, 3, 12246, 26100, 12576, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26610, 0, 3, 12252, 26110, 12594, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 26640, 0, 3, 12258, 26120, 12612, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 26670, 0, 3, 12270, 26130, 3462, 3480,
                                                 12666, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26730, 0, 3, 12288, 26160, 3480, 3498,
                                                 12702, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26790, 0, 3, 12306, 26190, 3498, 3516,
                                                 12738, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26850, 0, 3, 12324, 26220, 3516, 3534,
                                                 12774, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26910, 0, 3, 12342, 26250, 3534, 3552,
                                                 12810, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 26970, 0, 3, 12360, 26280, 3552, 3570,
                                                 12846, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27030, 0, 3, 12378, 26310, 3570, 3588,
                                                 12882, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27090, 0, 3, 12396, 26340, 3588, 3606,
                                                 12918, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27150, 0, 3, 12414, 26370, 3606, 3624,
                                                 12954, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27210, 0, 3, 12450, 26400, 3660, 3678,
                                                 13026, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27270, 0, 3, 12468, 26430, 3678, 3696,
                                                 13062, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27330, 0, 3, 12486, 26460, 3696, 3714,
                                                 13098, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27390, 0, 3, 12504, 26490, 3714, 3732,
                                                 13134, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27450, 0, 3, 12522, 26520, 3732, 3750,
                                                 13170, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27510, 0, 3, 12540, 26550, 3750, 3768,
                                                 13206, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27570, 0, 3, 12558, 26580, 3768, 3786,
                                                 13242, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27630, 0, 3, 12576, 26610, 3786, 3804,
                                                 13278, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 27690, 0, 3, 12594, 26640, 3804, 3822,
                                                 13314, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27750, 0, 3, 12630, 26670, 3858, 3888,
                                                 13350, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27850, 0, 3, 12666, 26730, 3888, 3918,
                                                 13410, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 27950, 0, 3, 12702, 26790, 3918, 3948,
                                                 13470, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28050, 0, 3, 12738, 26850, 3948, 3978,
                                                 13530, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28150, 0, 3, 12774, 26910, 3978, 4008,
                                                 13590, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28250, 0, 3, 12810, 26970, 4008, 4038,
                                                 13650, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28350, 0, 3, 12846, 27030, 4038, 4068,
                                                 13710, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28450, 0, 3, 12882, 27090, 4068, 4098,
                                                 13770, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28550, 0, 3, 12918, 27150, 4098, 4128,
                                                 13830, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28650, 0, 3, 12990, 27210, 4188, 4218,
                                                 13890, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28750, 0, 3, 13026, 27270, 4218, 4248,
                                                 13950, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28850, 0, 3, 13062, 27330, 4248, 4278,
                                                 14010, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 28950, 0, 3, 13098, 27390, 4278, 4308,
                                                 14070, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29050, 0, 3, 13134, 27450, 4308, 4338,
                                                 14130, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29150, 0, 3, 13170, 27510, 4338, 4368,
                                                 14190, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29250, 0, 3, 13206, 27570, 4368, 4398,
                                                 14250, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29350, 0, 3, 13242, 27630, 4398, 4428,
                                                 14310, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 29450, 0, 3, 13278, 27690, 4428, 4458,
                                                 14370, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29550, 0, 3, 26670, 26730, 13410, 27950,
                                                 4518, 4563, 14520, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29700, 0, 3, 26730, 26790, 13470, 28050,
                                                 4563, 4608, 14610, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 29850, 0, 3, 26790, 26850, 13530, 28150,
                                                 4608, 4653, 14700, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30000, 0, 3, 26850, 26910, 13590, 28250,
                                                 4653, 4698, 14790, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30150, 0, 3, 26910, 26970, 13650, 28350,
                                                 4698, 4743, 14880, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30300, 0, 3, 26970, 27030, 13710, 28450,
                                                 4743, 4788, 14970, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30450, 0, 3, 27030, 27090, 13770, 28550,
                                                 4788, 4833, 15060, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30600, 0, 3, 27210, 27270, 13950, 28850,
                                                 4923, 4968, 15240, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30750, 0, 3, 27270, 27330, 14010, 28950,
                                                 4968, 5013, 15330, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 30900, 0, 3, 27330, 27390, 14070, 29050,
                                                 5013, 5058, 15420, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31050, 0, 3, 27390, 27450, 14130, 29150,
                                                 5058, 5103, 15510, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31200, 0, 3, 27450, 27510, 14190, 29250,
                                                 5103, 5148, 15600, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31350, 0, 3, 27510, 27570, 14250, 29350,
                                                 5148, 5193, 15690, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 31500, 0, 3, 27570, 27630, 14310, 29450,
                                                 5193, 5238, 15780, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 31650, 0, 3, 27750, 27850, 14430, 29550,
                                                 5328, 5391, 15870, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 31860, 0, 3, 27850, 27950, 14520, 29700,
                                                 5391, 5454, 15996, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32070, 0, 3, 27950, 28050, 14610, 29850,
                                                 5454, 5517, 16122, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32280, 0, 3, 28050, 28150, 14700, 30000,
                                                 5517, 5580, 16248, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32490, 0, 3, 28150, 28250, 14790, 30150,
                                                 5580, 5643, 16374, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32700, 0, 3, 28250, 28350, 14880, 30300,
                                                 5643, 5706, 16500, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 32910, 0, 3, 28350, 28450, 14970, 30450,
                                                 5706, 5769, 16626, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33120, 0, 3, 28650, 28750, 15150, 30600,
                                                 5895, 5958, 16752, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33330, 0, 3, 28750, 28850, 15240, 30750,
                                                 5958, 6021, 16878, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33540, 0, 3, 28850, 28950, 15330, 30900,
                                                 6021, 6084, 17004, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33750, 0, 3, 28950, 29050, 15420, 31050,
                                                 6084, 6147, 17130, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 33960, 0, 3, 29050, 29150, 15510, 31200,
                                                 6147, 6210, 17256, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34170, 0, 3, 29150, 29250, 15600, 31350,
                                                 6210, 6273, 17382, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 34380, 0, 3, 29250, 29350, 15690, 31500,
                                                 6273, 6336, 17508, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 34590, 0, 3, 29550, 29700, 15996, 32070,
                                                 6462, 6546, 17802, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 34870, 0, 3, 29700, 29850, 16122, 32280,
                                                 6546, 6630, 17970, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35150, 0, 3, 29850, 30000, 16248, 32490,
                                                 6630, 6714, 18138, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35430, 0, 3, 30000, 30150, 16374, 32700,
                                                 6714, 6798, 18306, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35710, 0, 3, 30150, 30300, 16500, 32910,
                                                 6798, 6882, 18474, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 35990, 0, 3, 30600, 30750, 16878, 33540,
                                                 7050, 7134, 18810, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36270, 0, 3, 30750, 30900, 17004, 33750,
                                                 7134, 7218, 18978, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36550, 0, 3, 30900, 31050, 17130, 33960,
                                                 7218, 7302, 19146, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 36830, 0, 3, 31050, 31200, 17256, 34170,
                                                 7302, 7386, 19314, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 37110, 0, 3, 31200, 31350, 17382, 34380,
                                                 7386, 7470, 19482, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 37390, 0, 3, 31650, 31860, 17634, 34590,
                                                 7638, 7746, 19650, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 37750, 0, 3, 31860, 32070, 17802, 34870,
                                                 7746, 7854, 19866, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 38110, 0, 3, 32070, 32280, 17970, 35150,
                                                 7854, 7962, 20082, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 38470, 0, 3, 32280, 32490, 18138, 35430,
                                                 7962, 8070, 20298, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 38830, 0, 3, 32490, 32700, 18306, 35710,
                                                 8070, 8178, 20514, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39190, 0, 3, 33120, 33330, 18642, 35990,
                                                 8394, 8502, 20730, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39550, 0, 3, 33330, 33540, 18810, 36270,
                                                 8502, 8610, 20946, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 39910, 0, 3, 33540, 33750, 18978, 36550,
                                                 8610, 8718, 21162, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 40270, 0, 3, 33750, 33960, 19146, 36830,
                                                 8718, 8826, 21378, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 40630, 0, 3, 33960, 34170, 19314, 37110,
                                                 8826, 8934, 21594, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 40990, 0, 3, 34590, 34870, 19866, 38110,
                                                 9150, 9285, 22080, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 41440, 0, 3, 34870, 35150, 20082, 38470,
                                                 9285, 9420, 22350, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 41890, 0, 3, 35150, 35430, 20298, 38830,
                                                 9420, 9555, 22620, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 42340, 0, 3, 35990, 36270, 20946, 39910,
                                                 9825, 9960, 23160, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 42790, 0, 3, 36270, 36550, 21162, 40270,
                                                 9960, 10095, 23430, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 43240, 0, 3, 36550, 36830, 21378, 40630,
                                                 10095, 10230, 23700, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 43690, 0, 3, 37390, 37750, 21810, 40990,
                                                 10500, 10665, 23970, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 44240, 0, 3, 37750, 38110, 22080, 41440,
                                                 10665, 10830, 24300, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 44790, 0, 3, 38110, 38470, 22350, 41890,
                                                 10830, 10995, 24630, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 45340, 0, 3, 39190, 39550, 22890, 42340,
                                                 11325, 11490, 24960, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 45890, 0, 3, 39550, 39910, 23160, 42790,
                                                 11490, 11655, 25290, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 46440, 0, 3, 39910, 40270, 23430, 43240,
                                                 11655, 11820, 25620, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 46990, 3, 12150, 12156, 25960, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47005, 3, 12156, 12162, 25970, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47020, 3, 12162, 12168, 25980, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47035, 3, 12168, 12174, 25990, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47050, 3, 12174, 12180, 26000, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47065, 3, 12180, 12186, 26010, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47080, 3, 12186, 12192, 26020, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47095, 3, 12192, 12198, 26030, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47110, 3, 12210, 12216, 26050, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47125, 3, 12216, 12222, 26060, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47140, 3, 12222, 12228, 26070, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47155, 3, 12228, 12234, 26080, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47170, 3, 12234, 12240, 26090, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47185, 3, 12240, 12246, 26100, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47200, 3, 12246, 12252, 26110, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 47215, 3, 12252, 12258, 26120, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 47230, 0, 3, 25950, 46990, 26160, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47275, 0, 3, 25960, 47005, 26190, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47320, 0, 3, 25970, 47020, 26220, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47365, 0, 3, 25980, 47035, 26250, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47410, 0, 3, 25990, 47050, 26280, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47455, 0, 3, 26000, 47065, 26310, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47500, 0, 3, 26010, 47080, 26340, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47545, 0, 3, 26020, 47095, 26370, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47590, 0, 3, 26040, 47110, 26430, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47635, 0, 3, 26050, 47125, 26460, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47680, 0, 3, 26060, 47140, 26490, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47725, 0, 3, 26070, 47155, 26520, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47770, 0, 3, 26080, 47170, 26550, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47815, 0, 3, 26090, 47185, 26580, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47860, 0, 3, 26100, 47200, 26610, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 47905, 0, 3, 26110, 47215, 26640, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 47950, 0, 3, 26130, 47230, 12630, 12666,
                                                 26730, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48040, 0, 3, 26160, 47275, 12666, 12702,
                                                 26790, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48130, 0, 3, 26190, 47320, 12702, 12738,
                                                 26850, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48220, 0, 3, 26220, 47365, 12738, 12774,
                                                 26910, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48310, 0, 3, 26250, 47410, 12774, 12810,
                                                 26970, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48400, 0, 3, 26280, 47455, 12810, 12846,
                                                 27030, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48490, 0, 3, 26310, 47500, 12846, 12882,
                                                 27090, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48580, 0, 3, 26340, 47545, 12882, 12918,
                                                 27150, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48670, 0, 3, 26400, 47590, 12990, 13026,
                                                 27270, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48760, 0, 3, 26430, 47635, 13026, 13062,
                                                 27330, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48850, 0, 3, 26460, 47680, 13062, 13098,
                                                 27390, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 48940, 0, 3, 26490, 47725, 13098, 13134,
                                                 27450, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49030, 0, 3, 26520, 47770, 13134, 13170,
                                                 27510, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49120, 0, 3, 26550, 47815, 13170, 13206,
                                                 27570, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49210, 0, 3, 26580, 47860, 13206, 13242,
                                                 27630, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 49300, 0, 3, 26610, 47905, 13242, 13278,
                                                 27690, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49390, 0, 3, 26730, 48040, 13350, 13410,
                                                 27950, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49540, 0, 3, 26790, 48130, 13410, 13470,
                                                 28050, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49690, 0, 3, 26850, 48220, 13470, 13530,
                                                 28150, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49840, 0, 3, 26910, 48310, 13530, 13590,
                                                 28250, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 49990, 0, 3, 26970, 48400, 13590, 13650,
                                                 28350, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50140, 0, 3, 27030, 48490, 13650, 13710,
                                                 28450, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50290, 0, 3, 27090, 48580, 13710, 13770,
                                                 28550, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50440, 0, 3, 27270, 48760, 13890, 13950,
                                                 28850, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50590, 0, 3, 27330, 48850, 13950, 14010,
                                                 28950, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50740, 0, 3, 27390, 48940, 14010, 14070,
                                                 29050, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 50890, 0, 3, 27450, 49030, 14070, 14130,
                                                 29150, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51040, 0, 3, 27510, 49120, 14130, 14190,
                                                 29250, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51190, 0, 3, 27570, 49210, 14190, 14250,
                                                 29350, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 51340, 0, 3, 27630, 49300, 14250, 14310,
                                                 29450, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 51490, 0, 3, 47950, 48040, 27950, 49540,
                                                 14430, 14520, 29700, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 51715, 0, 3, 48040, 48130, 28050, 49690,
                                                 14520, 14610, 29850, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 51940, 0, 3, 48130, 48220, 28150, 49840,
                                                 14610, 14700, 30000, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52165, 0, 3, 48220, 48310, 28250, 49990,
                                                 14700, 14790, 30150, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52390, 0, 3, 48310, 48400, 28350, 50140,
                                                 14790, 14880, 30300, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52615, 0, 3, 48400, 48490, 28450, 50290,
                                                 14880, 14970, 30450, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 52840, 0, 3, 48670, 48760, 28850, 50590,
                                                 15150, 15240, 30750, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53065, 0, 3, 48760, 48850, 28950, 50740,
                                                 15240, 15330, 30900, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53290, 0, 3, 48850, 48940, 29050, 50890,
                                                 15330, 15420, 31050, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53515, 0, 3, 48940, 49030, 29150, 51040,
                                                 15420, 15510, 31200, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53740, 0, 3, 49030, 49120, 29250, 51190,
                                                 15510, 15600, 31350, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 53965, 0, 3, 49120, 49210, 29350, 51340,
                                                 15600, 15690, 31500, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 54190, 0, 3, 49390, 49540, 29700, 51715,
                                                 15870, 15996, 32070, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 54505, 0, 3, 49540, 49690, 29850, 51940,
                                                 15996, 16122, 32280, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 54820, 0, 3, 49690, 49840, 30000, 52165,
                                                 16122, 16248, 32490, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55135, 0, 3, 49840, 49990, 30150, 52390,
                                                 16248, 16374, 32700, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55450, 0, 3, 49990, 50140, 30300, 52615,
                                                 16374, 16500, 32910, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 55765, 0, 3, 50440, 50590, 30750, 53065,
                                                 16752, 16878, 33540, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56080, 0, 3, 50590, 50740, 30900, 53290,
                                                 16878, 17004, 33750, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56395, 0, 3, 50740, 50890, 31050, 53515,
                                                 17004, 17130, 33960, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 56710, 0, 3, 50890, 51040, 31200, 53740,
                                                 17130, 17256, 34170, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 57025, 0, 3, 51040, 51190, 31350, 53965,
                                                 17256, 17382, 34380, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 57340, 0, 3, 51490, 51715, 32070, 54505,
                                                 17634, 17802, 34870, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 57760, 0, 3, 51715, 51940, 32280, 54820,
                                                 17802, 17970, 35150, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 58180, 0, 3, 51940, 52165, 32490, 55135,
                                                 17970, 18138, 35430, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 58600, 0, 3, 52165, 52390, 32700, 55450,
                                                 18138, 18306, 35710, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59020, 0, 3, 52840, 53065, 33540, 56080,
                                                 18642, 18810, 36270, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59440, 0, 3, 53065, 53290, 33750, 56395,
                                                 18810, 18978, 36550, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 59860, 0, 3, 53290, 53515, 33960, 56710,
                                                 18978, 19146, 36830, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 60280, 0, 3, 53515, 53740, 34170, 57025,
                                                 19146, 19314, 37110, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 60700, 0, 3, 54190, 54505, 34870, 57760,
                                                 19650, 19866, 38110, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 61240, 0, 3, 54505, 54820, 35150, 58180,
                                                 19866, 20082, 38470, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 61780, 0, 3, 54820, 55135, 35430, 58600,
                                                 20082, 20298, 38830, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 62320, 0, 3, 55765, 56080, 36270, 59440,
                                                 20730, 20946, 39910, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 62860, 0, 3, 56080, 56395, 36550, 59860,
                                                 20946, 21162, 40270, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 63400, 0, 3, 56395, 56710, 36830, 60280,
                                                 21162, 21378, 40630, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 63940, 0, 3, 57340, 57760, 38110, 61240,
                                                 21810, 22080, 41440, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 64615, 0, 3, 57760, 58180, 38470, 61780,
                                                 22080, 22350, 41890, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 65290, 0, 3, 59020, 59440, 39910, 62860,
                                                 22890, 23160, 42790, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 65965, 0, 3, 59440, 59860, 40270, 63400,
                                                 23160, 23430, 43240, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 66640, 0, 3, 60700, 61240, 41440, 64615,
                                                 23970, 24300, 44790, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 67465, 0, 3, 62320, 62860, 42790, 65965,
                                                 24960, 25290, 46440, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68290, 3, 25950, 25960, 47005, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68311, 3, 25960, 25970, 47020, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68332, 3, 25970, 25980, 47035, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68353, 3, 25980, 25990, 47050, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68374, 3, 25990, 26000, 47065, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68395, 3, 26000, 26010, 47080, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68416, 3, 26010, 26020, 47095, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68437, 3, 26040, 26050, 47125, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68458, 3, 26050, 26060, 47140, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68479, 3, 26060, 26070, 47155, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68500, 3, 26070, 26080, 47170, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68521, 3, 26080, 26090, 47185, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68542, 3, 26090, 26100, 47200, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 68563, 3, 26100, 26110, 47215, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 68584, 0, 3, 46990, 68290, 47275, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68647, 0, 3, 47005, 68311, 47320, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68710, 0, 3, 47020, 68332, 47365, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68773, 0, 3, 47035, 68353, 47410, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68836, 0, 3, 47050, 68374, 47455, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68899, 0, 3, 47065, 68395, 47500, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 68962, 0, 3, 47080, 68416, 47545, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69025, 0, 3, 47110, 68437, 47635, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69088, 0, 3, 47125, 68458, 47680, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69151, 0, 3, 47140, 68479, 47725, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69214, 0, 3, 47155, 68500, 47770, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69277, 0, 3, 47170, 68521, 47815, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69340, 0, 3, 47185, 68542, 47860, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 69403, 0, 3, 47200, 68563, 47905, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 69466, 0, 3, 47230, 68584, 26670, 26730,
                                                 48040, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 69592, 0, 3, 47275, 68647, 26730, 26790,
                                                 48130, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 69718, 0, 3, 47320, 68710, 26790, 26850,
                                                 48220, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 69844, 0, 3, 47365, 68773, 26850, 26910,
                                                 48310, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 69970, 0, 3, 47410, 68836, 26910, 26970,
                                                 48400, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70096, 0, 3, 47455, 68899, 26970, 27030,
                                                 48490, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70222, 0, 3, 47500, 68962, 27030, 27090,
                                                 48580, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70348, 0, 3, 47590, 69025, 27210, 27270,
                                                 48760, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70474, 0, 3, 47635, 69088, 27270, 27330,
                                                 48850, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70600, 0, 3, 47680, 69151, 27330, 27390,
                                                 48940, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70726, 0, 3, 47725, 69214, 27390, 27450,
                                                 49030, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70852, 0, 3, 47770, 69277, 27450, 27510,
                                                 49120, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 70978, 0, 3, 47815, 69340, 27510, 27570,
                                                 49210, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 71104, 0, 3, 47860, 69403, 27570, 27630,
                                                 49300, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 71230, 0, 3, 47950, 69466, 27750, 27850,
                                                 49390, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 71440, 0, 3, 48040, 69592, 27850, 27950,
                                                 49540, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 71650, 0, 3, 48130, 69718, 27950, 28050,
                                                 49690, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 71860, 0, 3, 48220, 69844, 28050, 28150,
                                                 49840, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 72070, 0, 3, 48310, 69970, 28150, 28250,
                                                 49990, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 72280, 0, 3, 48400, 70096, 28250, 28350,
                                                 50140, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 72490, 0, 3, 48490, 70222, 28350, 28450,
                                                 50290, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 72700, 0, 3, 48670, 70348, 28650, 28750,
                                                 50440, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 72910, 0, 3, 48760, 70474, 28750, 28850,
                                                 50590, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 73120, 0, 3, 48850, 70600, 28850, 28950,
                                                 50740, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 73330, 0, 3, 48940, 70726, 28950, 29050,
                                                 50890, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 73540, 0, 3, 49030, 70852, 29050, 29150,
                                                 51040, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 73750, 0, 3, 49120, 70978, 29150, 29250,
                                                 51190, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 73960, 0, 3, 49210, 71104, 29250, 29350,
                                                 51340, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 74170, 0, 3, 69466, 69592, 49540, 71650,
                                                 29550, 29700, 51715, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 74485, 0, 3, 69592, 69718, 49690, 71860,
                                                 29700, 29850, 51940, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 74800, 0, 3, 69718, 69844, 49840, 72070,
                                                 29850, 30000, 52165, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 75115, 0, 3, 69844, 69970, 49990, 72280,
                                                 30000, 30150, 52390, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 75430, 0, 3, 69970, 70096, 50140, 72490,
                                                 30150, 30300, 52615, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 75745, 0, 3, 70348, 70474, 50590, 73120,
                                                 30600, 30750, 53065, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 76060, 0, 3, 70474, 70600, 50740, 73330,
                                                 30750, 30900, 53290, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 76375, 0, 3, 70600, 70726, 50890, 73540,
                                                 30900, 31050, 53515, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 76690, 0, 3, 70726, 70852, 51040, 73750,
                                                 31050, 31200, 53740, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 77005, 0, 3, 70852, 70978, 51190, 73960,
                                                 31200, 31350, 53965, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 77320, 0, 3, 71230, 71440, 51490, 74170,
                                                 31650, 31860, 54190, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 77761, 0, 3, 71440, 71650, 51715, 74485,
                                                 31860, 32070, 54505, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 78202, 0, 3, 71650, 71860, 51940, 74800,
                                                 32070, 32280, 54820, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 78643, 0, 3, 71860, 72070, 52165, 75115,
                                                 32280, 32490, 55135, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 79084, 0, 3, 72070, 72280, 52390, 75430,
                                                 32490, 32700, 55450, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 79525, 0, 3, 72700, 72910, 52840, 75745,
                                                 33120, 33330, 55765, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 79966, 0, 3, 72910, 73120, 53065, 76060,
                                                 33330, 33540, 56080, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 80407, 0, 3, 73120, 73330, 53290, 76375,
                                                 33540, 33750, 56395, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 80848, 0, 3, 73330, 73540, 53515, 76690,
                                                 33750, 33960, 56710, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 81289, 0, 3, 73540, 73750, 53740, 77005,
                                                 33960, 34170, 57025, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 81730, 0, 3, 74170, 74485, 54505, 78202,
                                                 34590, 34870, 57760, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 82318, 0, 3, 74485, 74800, 54820, 78643,
                                                 34870, 35150, 58180, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 82906, 0, 3, 74800, 75115, 55135, 79084,
                                                 35150, 35430, 58600, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 83494, 0, 3, 75745, 76060, 56080, 80407,
                                                 35990, 36270, 59440, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 84082, 0, 3, 76060, 76375, 56395, 80848,
                                                 36270, 36550, 59860, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 84670, 0, 3, 76375, 76690, 56710, 81289,
                                                 36550, 36830, 60280, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 85258, 0, 3, 77320, 77761, 57340, 81730,
                                                 37390, 37750, 60700, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 86014, 0, 3, 77761, 78202, 57760, 82318,
                                                 37750, 38110, 61240, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 86770, 0, 3, 78202, 78643, 58180, 82906,
                                                 38110, 38470, 61780, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 87526, 0, 3, 79525, 79966, 59020, 83494,
                                                 39190, 39550, 62320, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 88282, 0, 3, 79966, 80407, 59440, 84082,
                                                 39550, 39910, 62860, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 89038, 0, 3, 80407, 80848, 59860, 84670,
                                                 39910, 40270, 63400, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 89794, 0, 3, 81730, 82318, 61240, 86770,
                                                 40990, 41440, 64615, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 90739, 0, 3, 83494, 84082, 62860, 89038,
                                                 42340, 42790, 65965, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 91684, 0, 3, 85258, 86014, 63940, 89794,
                                                 43690, 44240, 66640, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 92839, 0, 3, 87526, 88282, 65290, 90739,
                                                 45340, 45890, 67465, ncols, alpha, beta, p);

            simdgeo::geom_l_x(buffer, 93994, 87526, 92839, 1, 21, ncols, alpha);

            simdgeo::geom_l_y(buffer, 94939, 87526, 92839, 1, 21, ncols, alpha);

            simdgeo::geom_l_z(buffer, 95884, 87526, 92839, 1, 21, ncols, alpha);

            simdgeo::geom_l_x(buffer, 96829, 85258, 91684, 1, 21, ncols, alpha);

            simdgeo::geom_l_y(buffer, 97774, 85258, 91684, 1, 21, ncols, alpha);

            simdgeo::geom_l_z(buffer, 98719, 85258, 91684, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 99664, 96829, 2835, ncols);

            simdfunc::contract_primitives(buffer, 102499, 93994, 2835, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 105334, 102499, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 105334, 11, nmax);

    simdtrf::transform_h_inner(buffer, 105334, 103444, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 187 * nvalues, nvalues, buffer, 105334, 11, nmax);

    simdtrf::transform_h_inner(buffer, 105334, 104389, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 374 * nvalues, nvalues, buffer, 105334, 11, nmax);

    simdtrf::transform_h_inner(buffer, 105334, 99664, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 561 * nvalues, nvalues, buffer, 105334, 11, nmax);

    simdtrf::transform_h_inner(buffer, 105334, 100609, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 748 * nvalues, nvalues, buffer, 105334, 11, nmax);

    simdtrf::transform_h_inner(buffer, 105334, 101554, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 935 * nvalues, nvalues, buffer, 105334, 11, nmax);
}

}  // namespace simdt2ceri
