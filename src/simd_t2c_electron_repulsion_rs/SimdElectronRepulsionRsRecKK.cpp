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


#include "SimdElectronRepulsionRsRecKK.hpp"

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
compute_rs_kk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_kk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 112130, 108998, 2592, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 2108, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2111, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2114, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2117, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2120, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2123, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2126, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2129, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2132, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2135, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2138, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2141, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2144, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2147, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2150, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2153, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2156, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2159, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2162, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2165, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2168, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2171, 3, 35, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2174, 3, 9, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2183, 3, 10, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2192, 3, 11, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2201, 3, 12, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2210, 3, 13, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2219, 3, 14, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2228, 3, 15, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2237, 3, 16, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2246, 3, 17, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2255, 3, 18, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2264, 3, 19, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2273, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2282, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2291, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2300, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2309, 3, 28, 93, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2318, 3, 29, 96, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2327, 3, 30, 99, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2336, 3, 31, 102, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2345, 3, 32, 105, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2354, 3, 33, 108, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2363, 3, 34, 111, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2372, 0, 3, 39, 2174, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2390, 0, 3, 42, 2183, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2408, 0, 3, 45, 2192, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2426, 0, 3, 48, 2201, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2444, 0, 3, 51, 2210, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2462, 0, 3, 54, 2219, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2480, 0, 3, 57, 2228, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2498, 0, 3, 60, 2237, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2516, 0, 3, 63, 2246, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2534, 0, 3, 66, 2255, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2552, 0, 3, 69, 2264, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2570, 0, 3, 78, 2273, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2588, 0, 3, 81, 2282, 198, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2606, 0, 3, 84, 2291, 204, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2624, 0, 3, 87, 2300, 210, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2642, 0, 3, 90, 2309, 216, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2660, 0, 3, 93, 2318, 222, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2678, 0, 3, 96, 2327, 228, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2696, 0, 3, 99, 2336, 234, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2714, 0, 3, 102, 2345, 240, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2732, 0, 3, 105, 2354, 246, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2750, 0, 3, 108, 2363, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2768, 0, 3, 114, 2372, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2798, 0, 3, 120, 2390, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2828, 0, 3, 126, 2408, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2858, 0, 3, 132, 2426, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2888, 0, 3, 138, 2444, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2918, 0, 3, 144, 2462, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2948, 0, 3, 150, 2480, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2978, 0, 3, 156, 2498, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3008, 0, 3, 162, 2516, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3038, 0, 3, 168, 2534, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3068, 0, 3, 174, 2552, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3098, 0, 3, 186, 2570, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3128, 0, 3, 192, 2588, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3158, 0, 3, 198, 2606, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3188, 0, 3, 204, 2624, 398, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3218, 0, 3, 210, 2642, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3248, 0, 3, 216, 2660, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3278, 0, 3, 222, 2678, 428, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3308, 0, 3, 228, 2696, 438, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3338, 0, 3, 234, 2714, 448, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3368, 0, 3, 240, 2732, 458, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3398, 0, 3, 246, 2750, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3428, 0, 3, 268, 2828, 493, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3473, 0, 3, 278, 2858, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3518, 0, 3, 288, 2888, 523, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3563, 0, 3, 298, 2918, 538, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3608, 0, 3, 308, 2948, 553, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3653, 0, 3, 318, 2978, 568, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3698, 0, 3, 328, 3008, 583, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3743, 0, 3, 338, 3038, 598, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3788, 0, 3, 348, 3068, 613, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3833, 0, 3, 378, 3158, 643, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3878, 0, 3, 388, 3188, 658, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3923, 0, 3, 398, 3218, 673, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3968, 0, 3, 408, 3248, 688, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4013, 0, 3, 418, 3278, 703, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4058, 0, 3, 428, 3308, 718, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4103, 0, 3, 438, 3338, 733, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4148, 0, 3, 448, 3368, 748, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4193, 0, 3, 458, 3398, 763, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4238, 0, 3, 478, 3428, 778, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4301, 0, 3, 493, 3473, 799, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4364, 0, 3, 508, 3518, 820, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4427, 0, 3, 523, 3563, 841, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4490, 0, 3, 538, 3608, 862, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4553, 0, 3, 553, 3653, 883, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4616, 0, 3, 568, 3698, 904, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4679, 0, 3, 583, 3743, 925, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4742, 0, 3, 598, 3788, 946, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4805, 0, 3, 628, 3833, 967, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4868, 0, 3, 643, 3878, 988, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4931, 0, 3, 658, 3923, 1009, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4994, 0, 3, 673, 3968, 1030, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5057, 0, 3, 688, 4013, 1051, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5120, 0, 3, 703, 4058, 1072, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5183, 0, 3, 718, 4103, 1093, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5246, 0, 3, 733, 4148, 1114, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5309, 0, 3, 748, 4193, 1135, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5372, 0, 3, 799, 4364, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5456, 0, 3, 820, 4427, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5540, 0, 3, 841, 4490, 1240, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5624, 0, 3, 862, 4553, 1268, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5708, 0, 3, 883, 4616, 1296, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5792, 0, 3, 904, 4679, 1324, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5876, 0, 3, 925, 4742, 1352, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5960, 0, 3, 988, 4931, 1408, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6044, 0, 3, 1009, 4994, 1436, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6128, 0, 3, 1030, 5057, 1464, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6212, 0, 3, 1051, 5120, 1492, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6296, 0, 3, 1072, 5183, 1520, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6380, 0, 3, 1093, 5246, 1548, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6464, 0, 3, 1114, 5309, 1576, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6548, 0, 3, 1156, 5372, 1604, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6656, 0, 3, 1184, 5456, 1640, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6764, 0, 3, 1212, 5540, 1676, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6872, 0, 3, 1240, 5624, 1712, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6980, 0, 3, 1268, 5708, 1748, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7088, 0, 3, 1296, 5792, 1784, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7196, 0, 3, 1324, 5876, 1820, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7304, 0, 3, 1380, 5960, 1856, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7412, 0, 3, 1408, 6044, 1892, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7520, 0, 3, 1436, 6128, 1928, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7628, 0, 3, 1464, 6212, 1964, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7736, 0, 3, 1492, 6296, 2000, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7844, 0, 3, 1520, 6380, 2036, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7952, 0, 3, 1548, 6464, 2072, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8060, 3, 9, 10, 2111, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8066, 3, 10, 11, 2114, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8072, 3, 11, 12, 2117, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8078, 3, 12, 13, 2120, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8084, 3, 13, 14, 2123, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8090, 3, 14, 15, 2126, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8096, 3, 15, 16, 2129, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8102, 3, 16, 17, 2132, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8108, 3, 17, 18, 2135, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8114, 3, 18, 19, 2138, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8120, 3, 24, 25, 2144, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8126, 3, 25, 26, 2147, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8132, 3, 26, 27, 2150, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8138, 3, 27, 28, 2153, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8144, 3, 28, 29, 2156, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8150, 3, 29, 30, 2159, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8156, 3, 30, 31, 2162, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8162, 3, 31, 32, 2165, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8168, 3, 32, 33, 2168, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8174, 3, 33, 34, 2171, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 8180, 0, 3, 2108, 8060, 2183, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8198, 0, 3, 2111, 8066, 2192, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8216, 0, 3, 2114, 8072, 2201, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8234, 0, 3, 2117, 8078, 2210, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8252, 0, 3, 2120, 8084, 2219, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8270, 0, 3, 2123, 8090, 2228, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8288, 0, 3, 2126, 8096, 2237, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8306, 0, 3, 2129, 8102, 2246, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8324, 0, 3, 2132, 8108, 2255, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8342, 0, 3, 2135, 8114, 2264, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8360, 0, 3, 2141, 8120, 2282, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8378, 0, 3, 2144, 8126, 2291, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8396, 0, 3, 2147, 8132, 2300, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8414, 0, 3, 2150, 8138, 2309, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8432, 0, 3, 2153, 8144, 2318, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8450, 0, 3, 2156, 8150, 2327, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8468, 0, 3, 2159, 8156, 2336, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8486, 0, 3, 2162, 8162, 2345, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8504, 0, 3, 2165, 8168, 2354, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8522, 0, 3, 2168, 8174, 2363, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 8540, 0, 3, 2174, 8180, 114, 120, 2390,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8576, 0, 3, 2183, 8198, 120, 126, 2408,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8612, 0, 3, 2192, 8216, 126, 132, 2426,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8648, 0, 3, 2201, 8234, 132, 138, 2444,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8684, 0, 3, 2210, 8252, 138, 144, 2462,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8720, 0, 3, 2219, 8270, 144, 150, 2480,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8756, 0, 3, 2228, 8288, 150, 156, 2498,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8792, 0, 3, 2237, 8306, 156, 162, 2516,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8828, 0, 3, 2246, 8324, 162, 168, 2534,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8864, 0, 3, 2255, 8342, 168, 174, 2552,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8900, 0, 3, 2273, 8360, 186, 192, 2588,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8936, 0, 3, 2282, 8378, 192, 198, 2606,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8972, 0, 3, 2291, 8396, 198, 204, 2624,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9008, 0, 3, 2300, 8414, 204, 210, 2642,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9044, 0, 3, 2309, 8432, 210, 216, 2660,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9080, 0, 3, 2318, 8450, 216, 222, 2678,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9116, 0, 3, 2327, 8468, 222, 228, 2696,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9152, 0, 3, 2336, 8486, 228, 234, 2714,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9188, 0, 3, 2345, 8504, 234, 240, 2732,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9224, 0, 3, 2354, 8522, 240, 246, 2750,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9260, 0, 3, 2390, 8576, 258, 268, 2828,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9320, 0, 3, 2408, 8612, 268, 278, 2858,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9380, 0, 3, 2426, 8648, 278, 288, 2888,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9440, 0, 3, 2444, 8684, 288, 298, 2918,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9500, 0, 3, 2462, 8720, 298, 308, 2948,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9560, 0, 3, 2480, 8756, 308, 318, 2978,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9620, 0, 3, 2498, 8792, 318, 328, 3008,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9680, 0, 3, 2516, 8828, 328, 338, 3038,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9740, 0, 3, 2534, 8864, 338, 348, 3068,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9800, 0, 3, 2588, 8936, 368, 378, 3158,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9860, 0, 3, 2606, 8972, 378, 388, 3188,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9920, 0, 3, 2624, 9008, 388, 398, 3218,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9980, 0, 3, 2642, 9044, 398, 408, 3248,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10040, 0, 3, 2660, 9080, 408, 418, 3278,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10100, 0, 3, 2678, 9116, 418, 428, 3308,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10160, 0, 3, 2696, 9152, 428, 438, 3338,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10220, 0, 3, 2714, 9188, 438, 448, 3368,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10280, 0, 3, 2732, 9224, 448, 458, 3398,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10340, 0, 3, 8540, 8576, 2828, 9320,
                                                 478, 493, 3473, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10430, 0, 3, 8576, 8612, 2858, 9380,
                                                 493, 508, 3518, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10520, 0, 3, 8612, 8648, 2888, 9440,
                                                 508, 523, 3563, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10610, 0, 3, 8648, 8684, 2918, 9500,
                                                 523, 538, 3608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10700, 0, 3, 8684, 8720, 2948, 9560,
                                                 538, 553, 3653, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10790, 0, 3, 8720, 8756, 2978, 9620,
                                                 553, 568, 3698, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10880, 0, 3, 8756, 8792, 3008, 9680,
                                                 568, 583, 3743, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10970, 0, 3, 8792, 8828, 3038, 9740,
                                                 583, 598, 3788, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11060, 0, 3, 8900, 8936, 3158, 9860,
                                                 628, 643, 3878, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11150, 0, 3, 8936, 8972, 3188, 9920,
                                                 643, 658, 3923, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11240, 0, 3, 8972, 9008, 3218, 9980,
                                                 658, 673, 3968, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11330, 0, 3, 9008, 9044, 3248, 10040,
                                                 673, 688, 4013, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11420, 0, 3, 9044, 9080, 3278, 10100,
                                                 688, 703, 4058, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11510, 0, 3, 9080, 9116, 3308, 10160,
                                                 703, 718, 4103, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11600, 0, 3, 9116, 9152, 3338, 10220,
                                                 718, 733, 4148, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11690, 0, 3, 9152, 9188, 3368, 10280,
                                                 733, 748, 4193, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11780, 0, 3, 9260, 9320, 3473, 10430,
                                                 778, 799, 4364, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11906, 0, 3, 9320, 9380, 3518, 10520,
                                                 799, 820, 4427, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12032, 0, 3, 9380, 9440, 3563, 10610,
                                                 820, 841, 4490, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12158, 0, 3, 9440, 9500, 3608, 10700,
                                                 841, 862, 4553, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12284, 0, 3, 9500, 9560, 3653, 10790,
                                                 862, 883, 4616, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12410, 0, 3, 9560, 9620, 3698, 10880,
                                                 883, 904, 4679, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12536, 0, 3, 9620, 9680, 3743, 10970,
                                                 904, 925, 4742, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12662, 0, 3, 9800, 9860, 3878, 11150,
                                                 967, 988, 4931, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12788, 0, 3, 9860, 9920, 3923, 11240,
                                                 988, 1009, 4994, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12914, 0, 3, 9920, 9980, 3968, 11330,
                                                 1009, 1030, 5057, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13040, 0, 3, 9980, 10040, 4013, 11420,
                                                 1030, 1051, 5120, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13166, 0, 3, 10040, 10100, 4058, 11510,
                                                 1051, 1072, 5183, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13292, 0, 3, 10100, 10160, 4103, 11600,
                                                 1072, 1093, 5246, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13418, 0, 3, 10160, 10220, 4148, 11690,
                                                 1093, 1114, 5309, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13544, 0, 3, 10340, 10430, 4364, 11906,
                                                 1156, 1184, 5456, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13712, 0, 3, 10430, 10520, 4427, 12032,
                                                 1184, 1212, 5540, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13880, 0, 3, 10520, 10610, 4490, 12158,
                                                 1212, 1240, 5624, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14048, 0, 3, 10610, 10700, 4553, 12284,
                                                 1240, 1268, 5708, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14216, 0, 3, 10700, 10790, 4616, 12410,
                                                 1268, 1296, 5792, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14384, 0, 3, 10790, 10880, 4679, 12536,
                                                 1296, 1324, 5876, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14552, 0, 3, 11060, 11150, 4931, 12788,
                                                 1380, 1408, 6044, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14720, 0, 3, 11150, 11240, 4994, 12914,
                                                 1408, 1436, 6128, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14888, 0, 3, 11240, 11330, 5057, 13040,
                                                 1436, 1464, 6212, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15056, 0, 3, 11330, 11420, 5120, 13166,
                                                 1464, 1492, 6296, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15224, 0, 3, 11420, 11510, 5183, 13292,
                                                 1492, 1520, 6380, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15392, 0, 3, 11510, 11600, 5246, 13418,
                                                 1520, 1548, 6464, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15560, 0, 3, 11780, 11906, 5456, 13712,
                                                 1604, 1640, 6764, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15776, 0, 3, 11906, 12032, 5540, 13880,
                                                 1640, 1676, 6872, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15992, 0, 3, 12032, 12158, 5624, 14048,
                                                 1676, 1712, 6980, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16208, 0, 3, 12158, 12284, 5708, 14216,
                                                 1712, 1748, 7088, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16424, 0, 3, 12284, 12410, 5792, 14384,
                                                 1748, 1784, 7196, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16640, 0, 3, 12662, 12788, 6044, 14720,
                                                 1856, 1892, 7520, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16856, 0, 3, 12788, 12914, 6128, 14888,
                                                 1892, 1928, 7628, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17072, 0, 3, 12914, 13040, 6212, 15056,
                                                 1928, 1964, 7736, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17288, 0, 3, 13040, 13166, 6296, 15224,
                                                 1964, 2000, 7844, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17504, 0, 3, 13166, 13292, 6380, 15392,
                                                 2000, 2036, 7952, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17720, 3, 2108, 2111, 8066, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17730, 3, 2111, 2114, 8072, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17740, 3, 2114, 2117, 8078, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17750, 3, 2117, 2120, 8084, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17760, 3, 2120, 2123, 8090, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17770, 3, 2123, 2126, 8096, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17780, 3, 2126, 2129, 8102, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17790, 3, 2129, 2132, 8108, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17800, 3, 2132, 2135, 8114, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17810, 3, 2141, 2144, 8126, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17820, 3, 2144, 2147, 8132, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17830, 3, 2147, 2150, 8138, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17840, 3, 2150, 2153, 8144, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17850, 3, 2153, 2156, 8150, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17860, 3, 2156, 2159, 8156, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17870, 3, 2159, 2162, 8162, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17880, 3, 2162, 2165, 8168, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17890, 3, 2165, 2168, 8174, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 17900, 0, 3, 8060, 17720, 8198, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17930, 0, 3, 8066, 17730, 8216, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17960, 0, 3, 8072, 17740, 8234, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17990, 0, 3, 8078, 17750, 8252, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18020, 0, 3, 8084, 17760, 8270, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18050, 0, 3, 8090, 17770, 8288, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18080, 0, 3, 8096, 17780, 8306, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18110, 0, 3, 8102, 17790, 8324, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18140, 0, 3, 8108, 17800, 8342, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18170, 0, 3, 8120, 17810, 8378, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18200, 0, 3, 8126, 17820, 8396, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18230, 0, 3, 8132, 17830, 8414, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18260, 0, 3, 8138, 17840, 8432, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18290, 0, 3, 8144, 17850, 8450, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18320, 0, 3, 8150, 17860, 8468, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18350, 0, 3, 8156, 17870, 8486, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18380, 0, 3, 8162, 17880, 8504, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18410, 0, 3, 8168, 17890, 8522, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 18440, 0, 3, 8180, 17900, 2372, 2390,
                                                 8576, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18500, 0, 3, 8198, 17930, 2390, 2408,
                                                 8612, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18560, 0, 3, 8216, 17960, 2408, 2426,
                                                 8648, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18620, 0, 3, 8234, 17990, 2426, 2444,
                                                 8684, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18680, 0, 3, 8252, 18020, 2444, 2462,
                                                 8720, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18740, 0, 3, 8270, 18050, 2462, 2480,
                                                 8756, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18800, 0, 3, 8288, 18080, 2480, 2498,
                                                 8792, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18860, 0, 3, 8306, 18110, 2498, 2516,
                                                 8828, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18920, 0, 3, 8324, 18140, 2516, 2534,
                                                 8864, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18980, 0, 3, 8360, 18170, 2570, 2588,
                                                 8936, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19040, 0, 3, 8378, 18200, 2588, 2606,
                                                 8972, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19100, 0, 3, 8396, 18230, 2606, 2624,
                                                 9008, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19160, 0, 3, 8414, 18260, 2624, 2642,
                                                 9044, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19220, 0, 3, 8432, 18290, 2642, 2660,
                                                 9080, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19280, 0, 3, 8450, 18320, 2660, 2678,
                                                 9116, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19340, 0, 3, 8468, 18350, 2678, 2696,
                                                 9152, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19400, 0, 3, 8486, 18380, 2696, 2714,
                                                 9188, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19460, 0, 3, 8504, 18410, 2714, 2732,
                                                 9224, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19520, 0, 3, 8540, 18440, 2768, 2798,
                                                 9260, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19620, 0, 3, 8576, 18500, 2798, 2828,
                                                 9320, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19720, 0, 3, 8612, 18560, 2828, 2858,
                                                 9380, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19820, 0, 3, 8648, 18620, 2858, 2888,
                                                 9440, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19920, 0, 3, 8684, 18680, 2888, 2918,
                                                 9500, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20020, 0, 3, 8720, 18740, 2918, 2948,
                                                 9560, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20120, 0, 3, 8756, 18800, 2948, 2978,
                                                 9620, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20220, 0, 3, 8792, 18860, 2978, 3008,
                                                 9680, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20320, 0, 3, 8828, 18920, 3008, 3038,
                                                 9740, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20420, 0, 3, 8900, 18980, 3098, 3128,
                                                 9800, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20520, 0, 3, 8936, 19040, 3128, 3158,
                                                 9860, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20620, 0, 3, 8972, 19100, 3158, 3188,
                                                 9920, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20720, 0, 3, 9008, 19160, 3188, 3218,
                                                 9980, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20820, 0, 3, 9044, 19220, 3218, 3248,
                                                 10040, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20920, 0, 3, 9080, 19280, 3248, 3278,
                                                 10100, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21020, 0, 3, 9116, 19340, 3278, 3308,
                                                 10160, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21120, 0, 3, 9152, 19400, 3308, 3338,
                                                 10220, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21220, 0, 3, 9188, 19460, 3338, 3368,
                                                 10280, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21320, 0, 3, 18440, 18500, 9320, 19720,
                                                 3428, 3473, 10430, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21470, 0, 3, 18500, 18560, 9380, 19820,
                                                 3473, 3518, 10520, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21620, 0, 3, 18560, 18620, 9440, 19920,
                                                 3518, 3563, 10610, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21770, 0, 3, 18620, 18680, 9500, 20020,
                                                 3563, 3608, 10700, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21920, 0, 3, 18680, 18740, 9560, 20120,
                                                 3608, 3653, 10790, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22070, 0, 3, 18740, 18800, 9620, 20220,
                                                 3653, 3698, 10880, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22220, 0, 3, 18800, 18860, 9680, 20320,
                                                 3698, 3743, 10970, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22370, 0, 3, 18980, 19040, 9860, 20620,
                                                 3833, 3878, 11150, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22520, 0, 3, 19040, 19100, 9920, 20720,
                                                 3878, 3923, 11240, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22670, 0, 3, 19100, 19160, 9980, 20820,
                                                 3923, 3968, 11330, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22820, 0, 3, 19160, 19220, 10040, 20920,
                                                 3968, 4013, 11420, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22970, 0, 3, 19220, 19280, 10100, 21020,
                                                 4013, 4058, 11510, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23120, 0, 3, 19280, 19340, 10160, 21120,
                                                 4058, 4103, 11600, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23270, 0, 3, 19340, 19400, 10220, 21220,
                                                 4103, 4148, 11690, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23420, 0, 3, 19520, 19620, 10340, 21320,
                                                 4238, 4301, 11780, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23630, 0, 3, 19620, 19720, 10430, 21470,
                                                 4301, 4364, 11906, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23840, 0, 3, 19720, 19820, 10520, 21620,
                                                 4364, 4427, 12032, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24050, 0, 3, 19820, 19920, 10610, 21770,
                                                 4427, 4490, 12158, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24260, 0, 3, 19920, 20020, 10700, 21920,
                                                 4490, 4553, 12284, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24470, 0, 3, 20020, 20120, 10790, 22070,
                                                 4553, 4616, 12410, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24680, 0, 3, 20120, 20220, 10880, 22220,
                                                 4616, 4679, 12536, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24890, 0, 3, 20420, 20520, 11060, 22370,
                                                 4805, 4868, 12662, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25100, 0, 3, 20520, 20620, 11150, 22520,
                                                 4868, 4931, 12788, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25310, 0, 3, 20620, 20720, 11240, 22670,
                                                 4931, 4994, 12914, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25520, 0, 3, 20720, 20820, 11330, 22820,
                                                 4994, 5057, 13040, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25730, 0, 3, 20820, 20920, 11420, 22970,
                                                 5057, 5120, 13166, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25940, 0, 3, 20920, 21020, 11510, 23120,
                                                 5120, 5183, 13292, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 26150, 0, 3, 21020, 21120, 11600, 23270,
                                                 5183, 5246, 13418, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26360, 0, 3, 21320, 21470, 11906, 23840,
                                                 5372, 5456, 13712, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26640, 0, 3, 21470, 21620, 12032, 24050,
                                                 5456, 5540, 13880, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26920, 0, 3, 21620, 21770, 12158, 24260,
                                                 5540, 5624, 14048, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27200, 0, 3, 21770, 21920, 12284, 24470,
                                                 5624, 5708, 14216, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27480, 0, 3, 21920, 22070, 12410, 24680,
                                                 5708, 5792, 14384, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27760, 0, 3, 22370, 22520, 12788, 25310,
                                                 5960, 6044, 14720, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28040, 0, 3, 22520, 22670, 12914, 25520,
                                                 6044, 6128, 14888, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28320, 0, 3, 22670, 22820, 13040, 25730,
                                                 6128, 6212, 15056, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28600, 0, 3, 22820, 22970, 13166, 25940,
                                                 6212, 6296, 15224, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28880, 0, 3, 22970, 23120, 13292, 26150,
                                                 6296, 6380, 15392, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29160, 0, 3, 23420, 23630, 13544, 26360,
                                                 6548, 6656, 15560, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29520, 0, 3, 23630, 23840, 13712, 26640,
                                                 6656, 6764, 15776, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29880, 0, 3, 23840, 24050, 13880, 26920,
                                                 6764, 6872, 15992, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30240, 0, 3, 24050, 24260, 14048, 27200,
                                                 6872, 6980, 16208, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30600, 0, 3, 24260, 24470, 14216, 27480,
                                                 6980, 7088, 16424, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30960, 0, 3, 24890, 25100, 14552, 27760,
                                                 7304, 7412, 16640, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 31320, 0, 3, 25100, 25310, 14720, 28040,
                                                 7412, 7520, 16856, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 31680, 0, 3, 25310, 25520, 14888, 28320,
                                                 7520, 7628, 17072, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 32040, 0, 3, 25520, 25730, 15056, 28600,
                                                 7628, 7736, 17288, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 32400, 0, 3, 25730, 25940, 15224, 28880,
                                                 7736, 7844, 17504, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32760, 3, 8060, 8066, 17730, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32775, 3, 8066, 8072, 17740, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32790, 3, 8072, 8078, 17750, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32805, 3, 8078, 8084, 17760, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32820, 3, 8084, 8090, 17770, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32835, 3, 8090, 8096, 17780, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32850, 3, 8096, 8102, 17790, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32865, 3, 8102, 8108, 17800, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32880, 3, 8120, 8126, 17820, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32895, 3, 8126, 8132, 17830, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32910, 3, 8132, 8138, 17840, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32925, 3, 8138, 8144, 17850, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32940, 3, 8144, 8150, 17860, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32955, 3, 8150, 8156, 17870, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32970, 3, 8156, 8162, 17880, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32985, 3, 8162, 8168, 17890, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 33000, 0, 3, 17720, 32760, 17930, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33045, 0, 3, 17730, 32775, 17960, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33090, 0, 3, 17740, 32790, 17990, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33135, 0, 3, 17750, 32805, 18020, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33180, 0, 3, 17760, 32820, 18050, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33225, 0, 3, 17770, 32835, 18080, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33270, 0, 3, 17780, 32850, 18110, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33315, 0, 3, 17790, 32865, 18140, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33360, 0, 3, 17810, 32880, 18200, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33405, 0, 3, 17820, 32895, 18230, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33450, 0, 3, 17830, 32910, 18260, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33495, 0, 3, 17840, 32925, 18290, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33540, 0, 3, 17850, 32940, 18320, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33585, 0, 3, 17860, 32955, 18350, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33630, 0, 3, 17870, 32970, 18380, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 33675, 0, 3, 17880, 32985, 18410, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 33720, 0, 3, 17900, 33000, 8540, 8576,
                                                 18500, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33810, 0, 3, 17930, 33045, 8576, 8612,
                                                 18560, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33900, 0, 3, 17960, 33090, 8612, 8648,
                                                 18620, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33990, 0, 3, 17990, 33135, 8648, 8684,
                                                 18680, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34080, 0, 3, 18020, 33180, 8684, 8720,
                                                 18740, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34170, 0, 3, 18050, 33225, 8720, 8756,
                                                 18800, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34260, 0, 3, 18080, 33270, 8756, 8792,
                                                 18860, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34350, 0, 3, 18110, 33315, 8792, 8828,
                                                 18920, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34440, 0, 3, 18170, 33360, 8900, 8936,
                                                 19040, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34530, 0, 3, 18200, 33405, 8936, 8972,
                                                 19100, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34620, 0, 3, 18230, 33450, 8972, 9008,
                                                 19160, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34710, 0, 3, 18260, 33495, 9008, 9044,
                                                 19220, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34800, 0, 3, 18290, 33540, 9044, 9080,
                                                 19280, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34890, 0, 3, 18320, 33585, 9080, 9116,
                                                 19340, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34980, 0, 3, 18350, 33630, 9116, 9152,
                                                 19400, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35070, 0, 3, 18380, 33675, 9152, 9188,
                                                 19460, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35160, 0, 3, 18500, 33810, 9260, 9320,
                                                 19720, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35310, 0, 3, 18560, 33900, 9320, 9380,
                                                 19820, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35460, 0, 3, 18620, 33990, 9380, 9440,
                                                 19920, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35610, 0, 3, 18680, 34080, 9440, 9500,
                                                 20020, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35760, 0, 3, 18740, 34170, 9500, 9560,
                                                 20120, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35910, 0, 3, 18800, 34260, 9560, 9620,
                                                 20220, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36060, 0, 3, 18860, 34350, 9620, 9680,
                                                 20320, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36210, 0, 3, 19040, 34530, 9800, 9860,
                                                 20620, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36360, 0, 3, 19100, 34620, 9860, 9920,
                                                 20720, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36510, 0, 3, 19160, 34710, 9920, 9980,
                                                 20820, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36660, 0, 3, 19220, 34800, 9980, 10040,
                                                 20920, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36810, 0, 3, 19280, 34890, 10040, 10100,
                                                 21020, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36960, 0, 3, 19340, 34980, 10100, 10160,
                                                 21120, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37110, 0, 3, 19400, 35070, 10160, 10220,
                                                 21220, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 37260, 0, 3, 33720, 33810, 19720, 35310,
                                                 10340, 10430, 21470, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 37485, 0, 3, 33810, 33900, 19820, 35460,
                                                 10430, 10520, 21620, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 37710, 0, 3, 33900, 33990, 19920, 35610,
                                                 10520, 10610, 21770, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 37935, 0, 3, 33990, 34080, 20020, 35760,
                                                 10610, 10700, 21920, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38160, 0, 3, 34080, 34170, 20120, 35910,
                                                 10700, 10790, 22070, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38385, 0, 3, 34170, 34260, 20220, 36060,
                                                 10790, 10880, 22220, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38610, 0, 3, 34440, 34530, 20620, 36360,
                                                 11060, 11150, 22520, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38835, 0, 3, 34530, 34620, 20720, 36510,
                                                 11150, 11240, 22670, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39060, 0, 3, 34620, 34710, 20820, 36660,
                                                 11240, 11330, 22820, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39285, 0, 3, 34710, 34800, 20920, 36810,
                                                 11330, 11420, 22970, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39510, 0, 3, 34800, 34890, 21020, 36960,
                                                 11420, 11510, 23120, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39735, 0, 3, 34890, 34980, 21120, 37110,
                                                 11510, 11600, 23270, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 39960, 0, 3, 35160, 35310, 21470, 37485,
                                                 11780, 11906, 23840, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40275, 0, 3, 35310, 35460, 21620, 37710,
                                                 11906, 12032, 24050, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40590, 0, 3, 35460, 35610, 21770, 37935,
                                                 12032, 12158, 24260, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40905, 0, 3, 35610, 35760, 21920, 38160,
                                                 12158, 12284, 24470, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41220, 0, 3, 35760, 35910, 22070, 38385,
                                                 12284, 12410, 24680, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41535, 0, 3, 36210, 36360, 22520, 38835,
                                                 12662, 12788, 25310, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41850, 0, 3, 36360, 36510, 22670, 39060,
                                                 12788, 12914, 25520, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42165, 0, 3, 36510, 36660, 22820, 39285,
                                                 12914, 13040, 25730, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42480, 0, 3, 36660, 36810, 22970, 39510,
                                                 13040, 13166, 25940, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42795, 0, 3, 36810, 36960, 23120, 39735,
                                                 13166, 13292, 26150, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 43110, 0, 3, 37260, 37485, 23840, 40275,
                                                 13544, 13712, 26640, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 43530, 0, 3, 37485, 37710, 24050, 40590,
                                                 13712, 13880, 26920, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 43950, 0, 3, 37710, 37935, 24260, 40905,
                                                 13880, 14048, 27200, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44370, 0, 3, 37935, 38160, 24470, 41220,
                                                 14048, 14216, 27480, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44790, 0, 3, 38610, 38835, 25310, 41850,
                                                 14552, 14720, 28040, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 45210, 0, 3, 38835, 39060, 25520, 42165,
                                                 14720, 14888, 28320, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 45630, 0, 3, 39060, 39285, 25730, 42480,
                                                 14888, 15056, 28600, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 46050, 0, 3, 39285, 39510, 25940, 42795,
                                                 15056, 15224, 28880, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 46470, 0, 3, 39960, 40275, 26640, 43530,
                                                 15560, 15776, 29880, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 47010, 0, 3, 40275, 40590, 26920, 43950,
                                                 15776, 15992, 30240, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 47550, 0, 3, 40590, 40905, 27200, 44370,
                                                 15992, 16208, 30600, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 48090, 0, 3, 41535, 41850, 28040, 45210,
                                                 16640, 16856, 31680, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 48630, 0, 3, 41850, 42165, 28320, 45630,
                                                 16856, 17072, 32040, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 49170, 0, 3, 42165, 42480, 28600, 46050,
                                                 17072, 17288, 32400, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49710, 3, 17720, 17730, 32775, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49731, 3, 17730, 17740, 32790, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49752, 3, 17740, 17750, 32805, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49773, 3, 17750, 17760, 32820, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49794, 3, 17760, 17770, 32835, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49815, 3, 17770, 17780, 32850, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49836, 3, 17780, 17790, 32865, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49857, 3, 17810, 17820, 32895, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49878, 3, 17820, 17830, 32910, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49899, 3, 17830, 17840, 32925, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49920, 3, 17840, 17850, 32940, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49941, 3, 17850, 17860, 32955, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49962, 3, 17860, 17870, 32970, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49983, 3, 17870, 17880, 32985, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 50004, 0, 3, 32760, 49710, 33045, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50067, 0, 3, 32775, 49731, 33090, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50130, 0, 3, 32790, 49752, 33135, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50193, 0, 3, 32805, 49773, 33180, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50256, 0, 3, 32820, 49794, 33225, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50319, 0, 3, 32835, 49815, 33270, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50382, 0, 3, 32850, 49836, 33315, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50445, 0, 3, 32880, 49857, 33405, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50508, 0, 3, 32895, 49878, 33450, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50571, 0, 3, 32910, 49899, 33495, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50634, 0, 3, 32925, 49920, 33540, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50697, 0, 3, 32940, 49941, 33585, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50760, 0, 3, 32955, 49962, 33630, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50823, 0, 3, 32970, 49983, 33675, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 50886, 0, 3, 33000, 50004, 18440, 18500,
                                                 33810, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51012, 0, 3, 33045, 50067, 18500, 18560,
                                                 33900, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51138, 0, 3, 33090, 50130, 18560, 18620,
                                                 33990, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51264, 0, 3, 33135, 50193, 18620, 18680,
                                                 34080, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51390, 0, 3, 33180, 50256, 18680, 18740,
                                                 34170, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51516, 0, 3, 33225, 50319, 18740, 18800,
                                                 34260, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51642, 0, 3, 33270, 50382, 18800, 18860,
                                                 34350, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51768, 0, 3, 33360, 50445, 18980, 19040,
                                                 34530, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51894, 0, 3, 33405, 50508, 19040, 19100,
                                                 34620, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 52020, 0, 3, 33450, 50571, 19100, 19160,
                                                 34710, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 52146, 0, 3, 33495, 50634, 19160, 19220,
                                                 34800, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 52272, 0, 3, 33540, 50697, 19220, 19280,
                                                 34890, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 52398, 0, 3, 33585, 50760, 19280, 19340,
                                                 34980, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 52524, 0, 3, 33630, 50823, 19340, 19400,
                                                 35070, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52650, 0, 3, 33720, 50886, 19520, 19620,
                                                 35160, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52860, 0, 3, 33810, 51012, 19620, 19720,
                                                 35310, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53070, 0, 3, 33900, 51138, 19720, 19820,
                                                 35460, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53280, 0, 3, 33990, 51264, 19820, 19920,
                                                 35610, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53490, 0, 3, 34080, 51390, 19920, 20020,
                                                 35760, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53700, 0, 3, 34170, 51516, 20020, 20120,
                                                 35910, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53910, 0, 3, 34260, 51642, 20120, 20220,
                                                 36060, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 54120, 0, 3, 34440, 51768, 20420, 20520,
                                                 36210, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 54330, 0, 3, 34530, 51894, 20520, 20620,
                                                 36360, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 54540, 0, 3, 34620, 52020, 20620, 20720,
                                                 36510, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 54750, 0, 3, 34710, 52146, 20720, 20820,
                                                 36660, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 54960, 0, 3, 34800, 52272, 20820, 20920,
                                                 36810, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 55170, 0, 3, 34890, 52398, 20920, 21020,
                                                 36960, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 55380, 0, 3, 34980, 52524, 21020, 21120,
                                                 37110, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55590, 0, 3, 50886, 51012, 35310, 53070,
                                                 21320, 21470, 37485, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55905, 0, 3, 51012, 51138, 35460, 53280,
                                                 21470, 21620, 37710, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 56220, 0, 3, 51138, 51264, 35610, 53490,
                                                 21620, 21770, 37935, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 56535, 0, 3, 51264, 51390, 35760, 53700,
                                                 21770, 21920, 38160, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 56850, 0, 3, 51390, 51516, 35910, 53910,
                                                 21920, 22070, 38385, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 57165, 0, 3, 51768, 51894, 36360, 54540,
                                                 22370, 22520, 38835, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 57480, 0, 3, 51894, 52020, 36510, 54750,
                                                 22520, 22670, 39060, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 57795, 0, 3, 52020, 52146, 36660, 54960,
                                                 22670, 22820, 39285, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 58110, 0, 3, 52146, 52272, 36810, 55170,
                                                 22820, 22970, 39510, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 58425, 0, 3, 52272, 52398, 36960, 55380,
                                                 22970, 23120, 39735, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 58740, 0, 3, 52650, 52860, 37260, 55590,
                                                 23420, 23630, 39960, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 59181, 0, 3, 52860, 53070, 37485, 55905,
                                                 23630, 23840, 40275, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 59622, 0, 3, 53070, 53280, 37710, 56220,
                                                 23840, 24050, 40590, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 60063, 0, 3, 53280, 53490, 37935, 56535,
                                                 24050, 24260, 40905, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 60504, 0, 3, 53490, 53700, 38160, 56850,
                                                 24260, 24470, 41220, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 60945, 0, 3, 54120, 54330, 38610, 57165,
                                                 24890, 25100, 41535, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 61386, 0, 3, 54330, 54540, 38835, 57480,
                                                 25100, 25310, 41850, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 61827, 0, 3, 54540, 54750, 39060, 57795,
                                                 25310, 25520, 42165, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 62268, 0, 3, 54750, 54960, 39285, 58110,
                                                 25520, 25730, 42480, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 62709, 0, 3, 54960, 55170, 39510, 58425,
                                                 25730, 25940, 42795, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 63150, 0, 3, 55590, 55905, 40275, 59622,
                                                 26360, 26640, 43530, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 63738, 0, 3, 55905, 56220, 40590, 60063,
                                                 26640, 26920, 43950, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 64326, 0, 3, 56220, 56535, 40905, 60504,
                                                 26920, 27200, 44370, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 64914, 0, 3, 57165, 57480, 41850, 61827,
                                                 27760, 28040, 45210, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 65502, 0, 3, 57480, 57795, 42165, 62268,
                                                 28040, 28320, 45630, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 66090, 0, 3, 57795, 58110, 42480, 62709,
                                                 28320, 28600, 46050, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 66678, 0, 3, 58740, 59181, 43110, 63150,
                                                 29160, 29520, 46470, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 67434, 0, 3, 59181, 59622, 43530, 63738,
                                                 29520, 29880, 47010, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 68190, 0, 3, 59622, 60063, 43950, 64326,
                                                 29880, 30240, 47550, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 68946, 0, 3, 60945, 61386, 44790, 64914,
                                                 30960, 31320, 48090, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 69702, 0, 3, 61386, 61827, 45210, 65502,
                                                 31320, 31680, 48630, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 70458, 0, 3, 61827, 62268, 45630, 66090,
                                                 31680, 32040, 49170, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71214, 3, 32760, 32775, 49731, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71242, 3, 32775, 32790, 49752, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71270, 3, 32790, 32805, 49773, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71298, 3, 32805, 32820, 49794, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71326, 3, 32820, 32835, 49815, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71354, 3, 32835, 32850, 49836, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71382, 3, 32880, 32895, 49878, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71410, 3, 32895, 32910, 49899, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71438, 3, 32910, 32925, 49920, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71466, 3, 32925, 32940, 49941, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71494, 3, 32940, 32955, 49962, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71522, 3, 32955, 32970, 49983, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 71550, 0, 3, 49710, 71214, 50067, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71634, 0, 3, 49731, 71242, 50130, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71718, 0, 3, 49752, 71270, 50193, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71802, 0, 3, 49773, 71298, 50256, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71886, 0, 3, 49794, 71326, 50319, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71970, 0, 3, 49815, 71354, 50382, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72054, 0, 3, 49857, 71382, 50508, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72138, 0, 3, 49878, 71410, 50571, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72222, 0, 3, 49899, 71438, 50634, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72306, 0, 3, 49920, 71466, 50697, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72390, 0, 3, 49941, 71494, 50760, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72474, 0, 3, 49962, 71522, 50823, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 72558, 0, 3, 50004, 71550, 33720, 33810,
                                                 51012, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 72726, 0, 3, 50067, 71634, 33810, 33900,
                                                 51138, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 72894, 0, 3, 50130, 71718, 33900, 33990,
                                                 51264, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73062, 0, 3, 50193, 71802, 33990, 34080,
                                                 51390, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73230, 0, 3, 50256, 71886, 34080, 34170,
                                                 51516, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73398, 0, 3, 50319, 71970, 34170, 34260,
                                                 51642, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73566, 0, 3, 50445, 72054, 34440, 34530,
                                                 51894, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73734, 0, 3, 50508, 72138, 34530, 34620,
                                                 52020, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73902, 0, 3, 50571, 72222, 34620, 34710,
                                                 52146, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 74070, 0, 3, 50634, 72306, 34710, 34800,
                                                 52272, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 74238, 0, 3, 50697, 72390, 34800, 34890,
                                                 52398, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 74406, 0, 3, 50760, 72474, 34890, 34980,
                                                 52524, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 74574, 0, 3, 51012, 72726, 35160, 35310,
                                                 53070, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 74854, 0, 3, 51138, 72894, 35310, 35460,
                                                 53280, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75134, 0, 3, 51264, 73062, 35460, 35610,
                                                 53490, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75414, 0, 3, 51390, 73230, 35610, 35760,
                                                 53700, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75694, 0, 3, 51516, 73398, 35760, 35910,
                                                 53910, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75974, 0, 3, 51894, 73734, 36210, 36360,
                                                 54540, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 76254, 0, 3, 52020, 73902, 36360, 36510,
                                                 54750, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 76534, 0, 3, 52146, 74070, 36510, 36660,
                                                 54960, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 76814, 0, 3, 52272, 74238, 36660, 36810,
                                                 55170, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 77094, 0, 3, 52398, 74406, 36810, 36960,
                                                 55380, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 77374, 0, 3, 72558, 72726, 53070, 74854,
                                                 37260, 37485, 55905, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 77794, 0, 3, 72726, 72894, 53280, 75134,
                                                 37485, 37710, 56220, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 78214, 0, 3, 72894, 73062, 53490, 75414,
                                                 37710, 37935, 56535, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 78634, 0, 3, 73062, 73230, 53700, 75694,
                                                 37935, 38160, 56850, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 79054, 0, 3, 73566, 73734, 54540, 76254,
                                                 38610, 38835, 57480, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 79474, 0, 3, 73734, 73902, 54750, 76534,
                                                 38835, 39060, 57795, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 79894, 0, 3, 73902, 74070, 54960, 76814,
                                                 39060, 39285, 58110, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 80314, 0, 3, 74070, 74238, 55170, 77094,
                                                 39285, 39510, 58425, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 80734, 0, 3, 74574, 74854, 55905, 77794,
                                                 39960, 40275, 59622, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 81322, 0, 3, 74854, 75134, 56220, 78214,
                                                 40275, 40590, 60063, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 81910, 0, 3, 75134, 75414, 56535, 78634,
                                                 40590, 40905, 60504, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 82498, 0, 3, 75974, 76254, 57480, 79474,
                                                 41535, 41850, 61827, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 83086, 0, 3, 76254, 76534, 57795, 79894,
                                                 41850, 42165, 62268, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 83674, 0, 3, 76534, 76814, 58110, 80314,
                                                 42165, 42480, 62709, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 84262, 0, 3, 77374, 77794, 59622, 81322,
                                                 43110, 43530, 63738, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 85046, 0, 3, 77794, 78214, 60063, 81910,
                                                 43530, 43950, 64326, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 85830, 0, 3, 79054, 79474, 61827, 83086,
                                                 44790, 45210, 65502, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 86614, 0, 3, 79474, 79894, 62268, 83674,
                                                 45210, 45630, 66090, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 87398, 0, 3, 80734, 81322, 63738, 85046,
                                                 46470, 47010, 68190, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 88406, 0, 3, 82498, 83086, 65502, 86614,
                                                 48090, 48630, 70458, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89414, 3, 49710, 49731, 71242, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89450, 3, 49731, 49752, 71270, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89486, 3, 49752, 49773, 71298, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89522, 3, 49773, 49794, 71326, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89558, 3, 49794, 49815, 71354, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89594, 3, 49857, 49878, 71410, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89630, 3, 49878, 49899, 71438, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89666, 3, 49899, 49920, 71466, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89702, 3, 49920, 49941, 71494, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 89738, 3, 49941, 49962, 71522, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 89774, 0, 3, 71214, 89414, 71634, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 89882, 0, 3, 71242, 89450, 71718, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 89990, 0, 3, 71270, 89486, 71802, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90098, 0, 3, 71298, 89522, 71886, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90206, 0, 3, 71326, 89558, 71970, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90314, 0, 3, 71382, 89594, 72138, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90422, 0, 3, 71410, 89630, 72222, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90530, 0, 3, 71438, 89666, 72306, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90638, 0, 3, 71466, 89702, 72390, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 90746, 0, 3, 71494, 89738, 72474, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 90854, 0, 3, 71550, 89774, 50886, 51012,
                                                 72726, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 91070, 0, 3, 71634, 89882, 51012, 51138,
                                                 72894, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 91286, 0, 3, 71718, 89990, 51138, 51264,
                                                 73062, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 91502, 0, 3, 71802, 90098, 51264, 51390,
                                                 73230, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 91718, 0, 3, 71886, 90206, 51390, 51516,
                                                 73398, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 91934, 0, 3, 72054, 90314, 51768, 51894,
                                                 73734, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92150, 0, 3, 72138, 90422, 51894, 52020,
                                                 73902, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92366, 0, 3, 72222, 90530, 52020, 52146,
                                                 74070, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92582, 0, 3, 72306, 90638, 52146, 52272,
                                                 74238, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92798, 0, 3, 72390, 90746, 52272, 52398,
                                                 74406, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 93014, 0, 3, 72558, 90854, 52650, 52860,
                                                 74574, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 93374, 0, 3, 72726, 91070, 52860, 53070,
                                                 74854, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 93734, 0, 3, 72894, 91286, 53070, 53280,
                                                 75134, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94094, 0, 3, 73062, 91502, 53280, 53490,
                                                 75414, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94454, 0, 3, 73230, 91718, 53490, 53700,
                                                 75694, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94814, 0, 3, 73566, 91934, 54120, 54330,
                                                 75974, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95174, 0, 3, 73734, 92150, 54330, 54540,
                                                 76254, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95534, 0, 3, 73902, 92366, 54540, 54750,
                                                 76534, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95894, 0, 3, 74070, 92582, 54750, 54960,
                                                 76814, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 96254, 0, 3, 74238, 92798, 54960, 55170,
                                                 77094, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 96614, 0, 3, 90854, 91070, 74854, 93734,
                                                 55590, 55905, 77794, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 97154, 0, 3, 91070, 91286, 75134, 94094,
                                                 55905, 56220, 78214, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 97694, 0, 3, 91286, 91502, 75414, 94454,
                                                 56220, 56535, 78634, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 98234, 0, 3, 91934, 92150, 76254, 95534,
                                                 57165, 57480, 79474, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 98774, 0, 3, 92150, 92366, 76534, 95894,
                                                 57480, 57795, 79894, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 99314, 0, 3, 92366, 92582, 76814, 96254,
                                                 57795, 58110, 80314, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 99854, 0, 3, 93014, 93374, 77374, 96614,
                                                 58740, 59181, 80734, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 100610, 0, 3, 93374, 93734, 77794,
                                                 97154, 59181, 59622, 81322, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 101366, 0, 3, 93734, 94094, 78214,
                                                 97694, 59622, 60063, 81910, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 102122, 0, 3, 94814, 95174, 79054,
                                                 98234, 60945, 61386, 82498, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 102878, 0, 3, 95174, 95534, 79474,
                                                 98774, 61386, 61827, 83086, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 103634, 0, 3, 95534, 95894, 79894,
                                                 99314, 61827, 62268, 83674, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 104390, 0, 3, 96614, 97154, 81322,
                                                 101366, 63150, 63738, 85046, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 105398, 0, 3, 98234, 98774, 83086,
                                                 103634, 64914, 65502, 86614, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 106406, 0, 3, 99854, 100610, 84262,
                                                 104390, 66678, 67434, 87398, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 107702, 0, 3, 102122, 102878, 85830,
                                                 105398, 68946, 69702, 88406, ncols, alpha, beta,
                                                 p);

            simdfunc::contract_primitives(buffer, 108998, 106406, 2592, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 111590, 110294, 36, 1, nmax);

    simdtrf::transform_k_outer_tri(values, nvalues, buffer, 111590, nmax);

    simdtrf::transform_k_inner(buffer, 111590, 108998, 36, 1, nmax);

    simdtrf::transform_k_outer_tri(values + 225 * nvalues, nvalues, buffer, 111590, nmax);
}

}  // namespace simdt2ceri
