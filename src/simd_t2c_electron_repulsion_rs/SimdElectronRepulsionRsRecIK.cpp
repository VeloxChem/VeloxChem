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


#include "SimdElectronRepulsionRsRecIK.hpp"

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
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_ik_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ik_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 74586, 72150, 2016, nvalues);

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
                                                9, 10, 11, 12, 13}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 20, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 100, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 103, 0, 33, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 7, 8, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 8, 9, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 9, 10, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 10, 11, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 11, 12, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 12, 13, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 13, 14, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 14, 15, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 15, 16, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 16, 17, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 17, 18, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 21, 22, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 22, 23, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 23, 24, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 24, 25, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 25, 26, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 26, 27, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 27, 28, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 28, 29, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 29, 30, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 30, 31, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 232, 0, 31, 32, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 34, 37, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 37, 40, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 40, 43, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 43, 46, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 46, 49, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 49, 52, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 52, 55, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 55, 58, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 58, 61, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 61, 64, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 70, 73, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 73, 76, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 76, 79, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 79, 82, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 82, 85, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 85, 88, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 398, 0, 88, 91, 214, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 408, 0, 91, 94, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 418, 0, 94, 97, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 428, 0, 97, 100, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 438, 0, 106, 112, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 453, 0, 112, 118, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 468, 0, 118, 124, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 483, 0, 124, 130, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 498, 0, 130, 136, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 513, 0, 136, 142, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 528, 0, 142, 148, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 543, 0, 148, 154, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 558, 0, 154, 160, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 573, 0, 172, 178, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 588, 0, 178, 184, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 603, 0, 184, 190, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 618, 0, 190, 196, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 633, 0, 196, 202, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 648, 0, 202, 208, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 663, 0, 208, 214, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 678, 0, 214, 220, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 693, 0, 220, 226, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 708, 0, 238, 248, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 729, 0, 248, 258, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 750, 0, 258, 268, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 771, 0, 268, 278, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 792, 0, 278, 288, 513, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 813, 0, 288, 298, 528, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 834, 0, 298, 308, 543, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 855, 0, 308, 318, 558, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 876, 0, 338, 348, 588, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 897, 0, 348, 358, 603, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 918, 0, 358, 368, 618, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 939, 0, 368, 378, 633, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 960, 0, 378, 388, 648, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 981, 0, 388, 398, 663, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1002, 0, 398, 408, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1023, 0, 408, 418, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1044, 0, 438, 453, 729, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1072, 0, 453, 468, 750, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1100, 0, 468, 483, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1128, 0, 483, 498, 792, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1156, 0, 498, 513, 813, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1184, 0, 513, 528, 834, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1212, 0, 528, 543, 855, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1240, 0, 573, 588, 897, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1268, 0, 588, 603, 918, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1296, 0, 603, 618, 939, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1324, 0, 618, 633, 960, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1352, 0, 633, 648, 981, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1380, 0, 648, 663, 1002, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1408, 0, 663, 678, 1023, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1436, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1439, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1442, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1445, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1448, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1451, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1454, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1457, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1460, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1463, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1466, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1469, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1472, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1475, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1478, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1481, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1484, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1487, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1490, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1493, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1496, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1499, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1502, 3, 8, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1511, 3, 9, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1520, 3, 10, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1529, 3, 11, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1538, 3, 12, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1547, 3, 13, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1556, 3, 14, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1565, 3, 15, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1574, 3, 16, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1583, 3, 17, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1592, 3, 18, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1601, 3, 22, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1610, 3, 23, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1619, 3, 24, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1628, 3, 25, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1637, 3, 26, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1646, 3, 27, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1655, 3, 28, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1664, 3, 29, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1673, 3, 30, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1682, 3, 31, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1691, 3, 32, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1700, 0, 3, 34, 1502, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1718, 0, 3, 37, 1511, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1736, 0, 3, 40, 1520, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1754, 0, 3, 43, 1529, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1772, 0, 3, 46, 1538, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1790, 0, 3, 49, 1547, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1808, 0, 3, 52, 1556, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1826, 0, 3, 55, 1565, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1844, 0, 3, 58, 1574, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1862, 0, 3, 61, 1583, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1880, 0, 3, 64, 1592, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1898, 0, 3, 70, 1601, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1916, 0, 3, 73, 1610, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1934, 0, 3, 76, 1619, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1952, 0, 3, 79, 1628, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1970, 0, 3, 82, 1637, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1988, 0, 3, 85, 1646, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2006, 0, 3, 88, 1655, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2024, 0, 3, 91, 1664, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2042, 0, 3, 94, 1673, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2060, 0, 3, 97, 1682, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2078, 0, 3, 100, 1691, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2096, 0, 3, 112, 1736, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2126, 0, 3, 118, 1754, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2156, 0, 3, 124, 1772, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2186, 0, 3, 130, 1790, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2216, 0, 3, 136, 1808, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2246, 0, 3, 142, 1826, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2276, 0, 3, 148, 1844, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2306, 0, 3, 154, 1862, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2336, 0, 3, 160, 1880, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2366, 0, 3, 178, 1934, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2396, 0, 3, 184, 1952, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2426, 0, 3, 190, 1970, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2456, 0, 3, 196, 1988, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2486, 0, 3, 202, 2006, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2516, 0, 3, 208, 2024, 398, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2546, 0, 3, 214, 2042, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2576, 0, 3, 220, 2060, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2606, 0, 3, 226, 2078, 428, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2636, 0, 3, 238, 2096, 438, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2681, 0, 3, 248, 2126, 453, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2726, 0, 3, 258, 2156, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2771, 0, 3, 268, 2186, 483, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2816, 0, 3, 278, 2216, 498, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2861, 0, 3, 288, 2246, 513, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2906, 0, 3, 298, 2276, 528, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2951, 0, 3, 308, 2306, 543, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2996, 0, 3, 318, 2336, 558, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3041, 0, 3, 338, 2366, 573, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3086, 0, 3, 348, 2396, 588, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3131, 0, 3, 358, 2426, 603, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3176, 0, 3, 368, 2456, 618, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3221, 0, 3, 378, 2486, 633, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3266, 0, 3, 388, 2516, 648, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3311, 0, 3, 398, 2546, 663, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3356, 0, 3, 408, 2576, 678, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3401, 0, 3, 418, 2606, 693, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3446, 0, 3, 453, 2726, 729, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3509, 0, 3, 468, 2771, 750, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3572, 0, 3, 483, 2816, 771, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3635, 0, 3, 498, 2861, 792, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3698, 0, 3, 513, 2906, 813, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3761, 0, 3, 528, 2951, 834, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3824, 0, 3, 543, 2996, 855, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3887, 0, 3, 588, 3131, 897, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3950, 0, 3, 603, 3176, 918, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4013, 0, 3, 618, 3221, 939, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4076, 0, 3, 633, 3266, 960, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4139, 0, 3, 648, 3311, 981, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4202, 0, 3, 663, 3356, 1002, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4265, 0, 3, 678, 3401, 1023, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4328, 0, 3, 708, 3446, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4412, 0, 3, 729, 3509, 1072, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4496, 0, 3, 750, 3572, 1100, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4580, 0, 3, 771, 3635, 1128, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4664, 0, 3, 792, 3698, 1156, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4748, 0, 3, 813, 3761, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4832, 0, 3, 834, 3824, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4916, 0, 3, 876, 3887, 1240, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5000, 0, 3, 897, 3950, 1268, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5084, 0, 3, 918, 4013, 1296, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5168, 0, 3, 939, 4076, 1324, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5252, 0, 3, 960, 4139, 1352, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5336, 0, 3, 981, 4202, 1380, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5420, 0, 3, 1002, 4265, 1408, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 5504, 3, 8, 9, 1439, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 5510, 3, 9, 10, 1442, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5516, 3, 10, 11, 1445, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5522, 3, 11, 12, 1448, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5528, 3, 12, 13, 1451, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5534, 3, 13, 14, 1454, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5540, 3, 14, 15, 1457, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5546, 3, 15, 16, 1460, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5552, 3, 16, 17, 1463, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5558, 3, 17, 18, 1466, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5564, 3, 22, 23, 1472, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5570, 3, 23, 24, 1475, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5576, 3, 24, 25, 1478, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5582, 3, 25, 26, 1481, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5588, 3, 26, 27, 1484, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5594, 3, 27, 28, 1487, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5600, 3, 28, 29, 1490, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5606, 3, 29, 30, 1493, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5612, 3, 30, 31, 1496, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5618, 3, 31, 32, 1499, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 5624, 0, 3, 1436, 5504, 1511, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5642, 0, 3, 1439, 5510, 1520, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5660, 0, 3, 1442, 5516, 1529, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5678, 0, 3, 1445, 5522, 1538, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5696, 0, 3, 1448, 5528, 1547, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5714, 0, 3, 1451, 5534, 1556, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5732, 0, 3, 1454, 5540, 1565, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5750, 0, 3, 1457, 5546, 1574, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5768, 0, 3, 1460, 5552, 1583, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5786, 0, 3, 1463, 5558, 1592, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5804, 0, 3, 1469, 5564, 1610, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5822, 0, 3, 1472, 5570, 1619, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5840, 0, 3, 1475, 5576, 1628, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5858, 0, 3, 1478, 5582, 1637, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5876, 0, 3, 1481, 5588, 1646, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5894, 0, 3, 1484, 5594, 1655, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5912, 0, 3, 1487, 5600, 1664, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5930, 0, 3, 1490, 5606, 1673, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5948, 0, 3, 1493, 5612, 1682, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5966, 0, 3, 1496, 5618, 1691, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 5984, 0, 3, 1511, 5642, 106, 112, 1736,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6020, 0, 3, 1520, 5660, 112, 118, 1754,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6056, 0, 3, 1529, 5678, 118, 124, 1772,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6092, 0, 3, 1538, 5696, 124, 130, 1790,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6128, 0, 3, 1547, 5714, 130, 136, 1808,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6164, 0, 3, 1556, 5732, 136, 142, 1826,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6200, 0, 3, 1565, 5750, 142, 148, 1844,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6236, 0, 3, 1574, 5768, 148, 154, 1862,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6272, 0, 3, 1583, 5786, 154, 160, 1880,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6308, 0, 3, 1610, 5822, 172, 178, 1934,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6344, 0, 3, 1619, 5840, 178, 184, 1952,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6380, 0, 3, 1628, 5858, 184, 190, 1970,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6416, 0, 3, 1637, 5876, 190, 196, 1988,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6452, 0, 3, 1646, 5894, 196, 202, 2006,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6488, 0, 3, 1655, 5912, 202, 208, 2024,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6524, 0, 3, 1664, 5930, 208, 214, 2042,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6560, 0, 3, 1673, 5948, 214, 220, 2060,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6596, 0, 3, 1682, 5966, 220, 226, 2078,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6632, 0, 3, 1736, 6020, 238, 248, 2126,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6692, 0, 3, 1754, 6056, 248, 258, 2156,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6752, 0, 3, 1772, 6092, 258, 268, 2186,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6812, 0, 3, 1790, 6128, 268, 278, 2216,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6872, 0, 3, 1808, 6164, 278, 288, 2246,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6932, 0, 3, 1826, 6200, 288, 298, 2276,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6992, 0, 3, 1844, 6236, 298, 308, 2306,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7052, 0, 3, 1862, 6272, 308, 318, 2336,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7112, 0, 3, 1934, 6344, 338, 348, 2396,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7172, 0, 3, 1952, 6380, 348, 358, 2426,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7232, 0, 3, 1970, 6416, 358, 368, 2456,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7292, 0, 3, 1988, 6452, 368, 378, 2486,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7352, 0, 3, 2006, 6488, 378, 388, 2516,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7412, 0, 3, 2024, 6524, 388, 398, 2546,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7472, 0, 3, 2042, 6560, 398, 408, 2576,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7532, 0, 3, 2060, 6596, 408, 418, 2606,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7592, 0, 3, 5984, 6020, 2126, 6692, 438,
                                                 453, 2726, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7682, 0, 3, 6020, 6056, 2156, 6752, 453,
                                                 468, 2771, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7772, 0, 3, 6056, 6092, 2186, 6812, 468,
                                                 483, 2816, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7862, 0, 3, 6092, 6128, 2216, 6872, 483,
                                                 498, 2861, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7952, 0, 3, 6128, 6164, 2246, 6932, 498,
                                                 513, 2906, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8042, 0, 3, 6164, 6200, 2276, 6992, 513,
                                                 528, 2951, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8132, 0, 3, 6200, 6236, 2306, 7052, 528,
                                                 543, 2996, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8222, 0, 3, 6308, 6344, 2396, 7172, 573,
                                                 588, 3131, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8312, 0, 3, 6344, 6380, 2426, 7232, 588,
                                                 603, 3176, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8402, 0, 3, 6380, 6416, 2456, 7292, 603,
                                                 618, 3221, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8492, 0, 3, 6416, 6452, 2486, 7352, 618,
                                                 633, 3266, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8582, 0, 3, 6452, 6488, 2516, 7412, 633,
                                                 648, 3311, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8672, 0, 3, 6488, 6524, 2546, 7472, 648,
                                                 663, 3356, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8762, 0, 3, 6524, 6560, 2576, 7532, 663,
                                                 678, 3401, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8852, 0, 3, 6632, 6692, 2726, 7682, 708,
                                                 729, 3509, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8978, 0, 3, 6692, 6752, 2771, 7772, 729,
                                                 750, 3572, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9104, 0, 3, 6752, 6812, 2816, 7862, 750,
                                                 771, 3635, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9230, 0, 3, 6812, 6872, 2861, 7952, 771,
                                                 792, 3698, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9356, 0, 3, 6872, 6932, 2906, 8042, 792,
                                                 813, 3761, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9482, 0, 3, 6932, 6992, 2951, 8132, 813,
                                                 834, 3824, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9608, 0, 3, 7112, 7172, 3131, 8312, 876,
                                                 897, 3950, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9734, 0, 3, 7172, 7232, 3176, 8402, 897,
                                                 918, 4013, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9860, 0, 3, 7232, 7292, 3221, 8492, 918,
                                                 939, 4076, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9986, 0, 3, 7292, 7352, 3266, 8582, 939,
                                                 960, 4139, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10112, 0, 3, 7352, 7412, 3311, 8672,
                                                 960, 981, 4202, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10238, 0, 3, 7412, 7472, 3356, 8762,
                                                 981, 1002, 4265, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10364, 0, 3, 7592, 7682, 3509, 8978,
                                                 1044, 1072, 4496, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10532, 0, 3, 7682, 7772, 3572, 9104,
                                                 1072, 1100, 4580, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10700, 0, 3, 7772, 7862, 3635, 9230,
                                                 1100, 1128, 4664, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10868, 0, 3, 7862, 7952, 3698, 9356,
                                                 1128, 1156, 4748, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11036, 0, 3, 7952, 8042, 3761, 9482,
                                                 1156, 1184, 4832, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11204, 0, 3, 8222, 8312, 3950, 9734,
                                                 1240, 1268, 5084, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11372, 0, 3, 8312, 8402, 4013, 9860,
                                                 1268, 1296, 5168, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11540, 0, 3, 8402, 8492, 4076, 9986,
                                                 1296, 1324, 5252, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11708, 0, 3, 8492, 8582, 4139, 10112,
                                                 1324, 1352, 5336, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11876, 0, 3, 8582, 8672, 4202, 10238,
                                                 1352, 1380, 5420, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12044, 3, 1436, 1439, 5510, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12054, 3, 1439, 1442, 5516, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12064, 3, 1442, 1445, 5522, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12074, 3, 1445, 1448, 5528, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12084, 3, 1448, 1451, 5534, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12094, 3, 1451, 1454, 5540, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12104, 3, 1454, 1457, 5546, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12114, 3, 1457, 1460, 5552, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12124, 3, 1460, 1463, 5558, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12134, 3, 1469, 1472, 5570, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12144, 3, 1472, 1475, 5576, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12154, 3, 1475, 1478, 5582, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12164, 3, 1478, 1481, 5588, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12174, 3, 1481, 1484, 5594, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12184, 3, 1484, 1487, 5600, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12194, 3, 1487, 1490, 5606, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12204, 3, 1490, 1493, 5612, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12214, 3, 1493, 1496, 5618, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 12224, 0, 3, 5504, 12044, 5642, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12254, 0, 3, 5510, 12054, 5660, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12284, 0, 3, 5516, 12064, 5678, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12314, 0, 3, 5522, 12074, 5696, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12344, 0, 3, 5528, 12084, 5714, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12374, 0, 3, 5534, 12094, 5732, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12404, 0, 3, 5540, 12104, 5750, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12434, 0, 3, 5546, 12114, 5768, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12464, 0, 3, 5552, 12124, 5786, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12494, 0, 3, 5564, 12134, 5822, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12524, 0, 3, 5570, 12144, 5840, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12554, 0, 3, 5576, 12154, 5858, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12584, 0, 3, 5582, 12164, 5876, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12614, 0, 3, 5588, 12174, 5894, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12644, 0, 3, 5594, 12184, 5912, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12674, 0, 3, 5600, 12194, 5930, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12704, 0, 3, 5606, 12204, 5948, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12734, 0, 3, 5612, 12214, 5966, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 12764, 0, 3, 5624, 12224, 1700, 1718,
                                                 5984, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12824, 0, 3, 5642, 12254, 1718, 1736,
                                                 6020, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12884, 0, 3, 5660, 12284, 1736, 1754,
                                                 6056, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12944, 0, 3, 5678, 12314, 1754, 1772,
                                                 6092, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13004, 0, 3, 5696, 12344, 1772, 1790,
                                                 6128, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13064, 0, 3, 5714, 12374, 1790, 1808,
                                                 6164, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13124, 0, 3, 5732, 12404, 1808, 1826,
                                                 6200, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13184, 0, 3, 5750, 12434, 1826, 1844,
                                                 6236, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13244, 0, 3, 5768, 12464, 1844, 1862,
                                                 6272, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13304, 0, 3, 5804, 12494, 1898, 1916,
                                                 6308, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13364, 0, 3, 5822, 12524, 1916, 1934,
                                                 6344, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13424, 0, 3, 5840, 12554, 1934, 1952,
                                                 6380, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13484, 0, 3, 5858, 12584, 1952, 1970,
                                                 6416, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13544, 0, 3, 5876, 12614, 1970, 1988,
                                                 6452, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13604, 0, 3, 5894, 12644, 1988, 2006,
                                                 6488, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13664, 0, 3, 5912, 12674, 2006, 2024,
                                                 6524, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13724, 0, 3, 5930, 12704, 2024, 2042,
                                                 6560, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13784, 0, 3, 5948, 12734, 2042, 2060,
                                                 6596, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13844, 0, 3, 6020, 12884, 2096, 2126,
                                                 6692, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13944, 0, 3, 6056, 12944, 2126, 2156,
                                                 6752, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14044, 0, 3, 6092, 13004, 2156, 2186,
                                                 6812, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14144, 0, 3, 6128, 13064, 2186, 2216,
                                                 6872, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14244, 0, 3, 6164, 13124, 2216, 2246,
                                                 6932, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14344, 0, 3, 6200, 13184, 2246, 2276,
                                                 6992, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14444, 0, 3, 6236, 13244, 2276, 2306,
                                                 7052, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14544, 0, 3, 6344, 13424, 2366, 2396,
                                                 7172, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14644, 0, 3, 6380, 13484, 2396, 2426,
                                                 7232, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14744, 0, 3, 6416, 13544, 2426, 2456,
                                                 7292, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14844, 0, 3, 6452, 13604, 2456, 2486,
                                                 7352, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14944, 0, 3, 6488, 13664, 2486, 2516,
                                                 7412, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15044, 0, 3, 6524, 13724, 2516, 2546,
                                                 7472, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15144, 0, 3, 6560, 13784, 2546, 2576,
                                                 7532, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15244, 0, 3, 12764, 12824, 6632, 13844,
                                                 2636, 2681, 7592, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15394, 0, 3, 12824, 12884, 6692, 13944,
                                                 2681, 2726, 7682, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15544, 0, 3, 12884, 12944, 6752, 14044,
                                                 2726, 2771, 7772, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15694, 0, 3, 12944, 13004, 6812, 14144,
                                                 2771, 2816, 7862, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15844, 0, 3, 13004, 13064, 6872, 14244,
                                                 2816, 2861, 7952, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15994, 0, 3, 13064, 13124, 6932, 14344,
                                                 2861, 2906, 8042, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16144, 0, 3, 13124, 13184, 6992, 14444,
                                                 2906, 2951, 8132, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16294, 0, 3, 13304, 13364, 7112, 14544,
                                                 3041, 3086, 8222, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16444, 0, 3, 13364, 13424, 7172, 14644,
                                                 3086, 3131, 8312, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16594, 0, 3, 13424, 13484, 7232, 14744,
                                                 3131, 3176, 8402, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16744, 0, 3, 13484, 13544, 7292, 14844,
                                                 3176, 3221, 8492, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16894, 0, 3, 13544, 13604, 7352, 14944,
                                                 3221, 3266, 8582, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17044, 0, 3, 13604, 13664, 7412, 15044,
                                                 3266, 3311, 8672, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17194, 0, 3, 13664, 13724, 7472, 15144,
                                                 3311, 3356, 8762, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17344, 0, 3, 13844, 13944, 7682, 15544,
                                                 3446, 3509, 8978, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17554, 0, 3, 13944, 14044, 7772, 15694,
                                                 3509, 3572, 9104, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17764, 0, 3, 14044, 14144, 7862, 15844,
                                                 3572, 3635, 9230, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17974, 0, 3, 14144, 14244, 7952, 15994,
                                                 3635, 3698, 9356, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18184, 0, 3, 14244, 14344, 8042, 16144,
                                                 3698, 3761, 9482, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18394, 0, 3, 14544, 14644, 8312, 16594,
                                                 3887, 3950, 9734, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18604, 0, 3, 14644, 14744, 8402, 16744,
                                                 3950, 4013, 9860, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18814, 0, 3, 14744, 14844, 8492, 16894,
                                                 4013, 4076, 9986, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19024, 0, 3, 14844, 14944, 8582, 17044,
                                                 4076, 4139, 10112, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19234, 0, 3, 14944, 15044, 8672, 17194,
                                                 4139, 4202, 10238, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19444, 0, 3, 15244, 15394, 8852, 17344,
                                                 4328, 4412, 10364, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19724, 0, 3, 15394, 15544, 8978, 17554,
                                                 4412, 4496, 10532, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20004, 0, 3, 15544, 15694, 9104, 17764,
                                                 4496, 4580, 10700, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20284, 0, 3, 15694, 15844, 9230, 17974,
                                                 4580, 4664, 10868, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20564, 0, 3, 15844, 15994, 9356, 18184,
                                                 4664, 4748, 11036, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20844, 0, 3, 16294, 16444, 9608, 18394,
                                                 4916, 5000, 11204, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21124, 0, 3, 16444, 16594, 9734, 18604,
                                                 5000, 5084, 11372, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21404, 0, 3, 16594, 16744, 9860, 18814,
                                                 5084, 5168, 11540, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21684, 0, 3, 16744, 16894, 9986, 19024,
                                                 5168, 5252, 11708, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21964, 0, 3, 16894, 17044, 10112, 19234,
                                                 5252, 5336, 11876, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22244, 3, 5504, 5510, 12054, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22259, 3, 5510, 5516, 12064, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22274, 3, 5516, 5522, 12074, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22289, 3, 5522, 5528, 12084, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22304, 3, 5528, 5534, 12094, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22319, 3, 5534, 5540, 12104, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22334, 3, 5540, 5546, 12114, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22349, 3, 5546, 5552, 12124, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22364, 3, 5564, 5570, 12144, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22379, 3, 5570, 5576, 12154, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22394, 3, 5576, 5582, 12164, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22409, 3, 5582, 5588, 12174, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22424, 3, 5588, 5594, 12184, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22439, 3, 5594, 5600, 12194, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22454, 3, 5600, 5606, 12204, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 22469, 3, 5606, 5612, 12214, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 22484, 0, 3, 12044, 22244, 12254, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22529, 0, 3, 12054, 22259, 12284, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22574, 0, 3, 12064, 22274, 12314, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22619, 0, 3, 12074, 22289, 12344, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22664, 0, 3, 12084, 22304, 12374, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22709, 0, 3, 12094, 22319, 12404, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22754, 0, 3, 12104, 22334, 12434, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22799, 0, 3, 12114, 22349, 12464, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22844, 0, 3, 12134, 22364, 12524, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22889, 0, 3, 12144, 22379, 12554, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22934, 0, 3, 12154, 22394, 12584, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 22979, 0, 3, 12164, 22409, 12614, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23024, 0, 3, 12174, 22424, 12644, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23069, 0, 3, 12184, 22439, 12674, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23114, 0, 3, 12194, 22454, 12704, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23159, 0, 3, 12204, 22469, 12734, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 23204, 0, 3, 12254, 22529, 5984, 6020,
                                                 12884, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23294, 0, 3, 12284, 22574, 6020, 6056,
                                                 12944, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23384, 0, 3, 12314, 22619, 6056, 6092,
                                                 13004, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23474, 0, 3, 12344, 22664, 6092, 6128,
                                                 13064, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23564, 0, 3, 12374, 22709, 6128, 6164,
                                                 13124, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23654, 0, 3, 12404, 22754, 6164, 6200,
                                                 13184, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23744, 0, 3, 12434, 22799, 6200, 6236,
                                                 13244, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23834, 0, 3, 12524, 22889, 6308, 6344,
                                                 13424, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23924, 0, 3, 12554, 22934, 6344, 6380,
                                                 13484, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24014, 0, 3, 12584, 22979, 6380, 6416,
                                                 13544, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24104, 0, 3, 12614, 23024, 6416, 6452,
                                                 13604, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24194, 0, 3, 12644, 23069, 6452, 6488,
                                                 13664, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24284, 0, 3, 12674, 23114, 6488, 6524,
                                                 13724, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24374, 0, 3, 12704, 23159, 6524, 6560,
                                                 13784, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24464, 0, 3, 12884, 23294, 6632, 6692,
                                                 13944, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24614, 0, 3, 12944, 23384, 6692, 6752,
                                                 14044, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24764, 0, 3, 13004, 23474, 6752, 6812,
                                                 14144, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24914, 0, 3, 13064, 23564, 6812, 6872,
                                                 14244, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25064, 0, 3, 13124, 23654, 6872, 6932,
                                                 14344, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25214, 0, 3, 13184, 23744, 6932, 6992,
                                                 14444, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25364, 0, 3, 13424, 23924, 7112, 7172,
                                                 14644, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25514, 0, 3, 13484, 24014, 7172, 7232,
                                                 14744, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25664, 0, 3, 13544, 24104, 7232, 7292,
                                                 14844, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25814, 0, 3, 13604, 24194, 7292, 7352,
                                                 14944, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25964, 0, 3, 13664, 24284, 7352, 7412,
                                                 15044, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 26114, 0, 3, 13724, 24374, 7412, 7472,
                                                 15144, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26264, 0, 3, 23204, 23294, 13944, 24614,
                                                 7592, 7682, 15544, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26489, 0, 3, 23294, 23384, 14044, 24764,
                                                 7682, 7772, 15694, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26714, 0, 3, 23384, 23474, 14144, 24914,
                                                 7772, 7862, 15844, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26939, 0, 3, 23474, 23564, 14244, 25064,
                                                 7862, 7952, 15994, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27164, 0, 3, 23564, 23654, 14344, 25214,
                                                 7952, 8042, 16144, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27389, 0, 3, 23834, 23924, 14644, 25514,
                                                 8222, 8312, 16594, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27614, 0, 3, 23924, 24014, 14744, 25664,
                                                 8312, 8402, 16744, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27839, 0, 3, 24014, 24104, 14844, 25814,
                                                 8402, 8492, 16894, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 28064, 0, 3, 24104, 24194, 14944, 25964,
                                                 8492, 8582, 17044, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 28289, 0, 3, 24194, 24284, 15044, 26114,
                                                 8582, 8672, 17194, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28514, 0, 3, 24464, 24614, 15544, 26489,
                                                 8852, 8978, 17554, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28829, 0, 3, 24614, 24764, 15694, 26714,
                                                 8978, 9104, 17764, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29144, 0, 3, 24764, 24914, 15844, 26939,
                                                 9104, 9230, 17974, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29459, 0, 3, 24914, 25064, 15994, 27164,
                                                 9230, 9356, 18184, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29774, 0, 3, 25364, 25514, 16594, 27614,
                                                 9608, 9734, 18604, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 30089, 0, 3, 25514, 25664, 16744, 27839,
                                                 9734, 9860, 18814, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 30404, 0, 3, 25664, 25814, 16894, 28064,
                                                 9860, 9986, 19024, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 30719, 0, 3, 25814, 25964, 17044, 28289,
                                                 9986, 10112, 19234, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 31034, 0, 3, 26264, 26489, 17554, 28829,
                                                 10364, 10532, 20004, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 31454, 0, 3, 26489, 26714, 17764, 29144,
                                                 10532, 10700, 20284, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 31874, 0, 3, 26714, 26939, 17974, 29459,
                                                 10700, 10868, 20564, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 32294, 0, 3, 27389, 27614, 18604, 30089,
                                                 11204, 11372, 21404, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 32714, 0, 3, 27614, 27839, 18814, 30404,
                                                 11372, 11540, 21684, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 33134, 0, 3, 27839, 28064, 19024, 30719,
                                                 11540, 11708, 21964, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33554, 3, 12044, 12054, 22259, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33575, 3, 12054, 12064, 22274, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33596, 3, 12064, 12074, 22289, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33617, 3, 12074, 12084, 22304, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33638, 3, 12084, 12094, 22319, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33659, 3, 12094, 12104, 22334, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33680, 3, 12104, 12114, 22349, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33701, 3, 12134, 12144, 22379, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33722, 3, 12144, 12154, 22394, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33743, 3, 12154, 12164, 22409, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33764, 3, 12164, 12174, 22424, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33785, 3, 12174, 12184, 22439, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33806, 3, 12184, 12194, 22454, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33827, 3, 12194, 12204, 22469, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 33848, 0, 3, 22244, 33554, 22529, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33911, 0, 3, 22259, 33575, 22574, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33974, 0, 3, 22274, 33596, 22619, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34037, 0, 3, 22289, 33617, 22664, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34100, 0, 3, 22304, 33638, 22709, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34163, 0, 3, 22319, 33659, 22754, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34226, 0, 3, 22334, 33680, 22799, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34289, 0, 3, 22364, 33701, 22889, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34352, 0, 3, 22379, 33722, 22934, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34415, 0, 3, 22394, 33743, 22979, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34478, 0, 3, 22409, 33764, 23024, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34541, 0, 3, 22424, 33785, 23069, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34604, 0, 3, 22439, 33806, 23114, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34667, 0, 3, 22454, 33827, 23159, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 34730, 0, 3, 22484, 33848, 12764, 12824,
                                                 23204, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34856, 0, 3, 22529, 33911, 12824, 12884,
                                                 23294, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34982, 0, 3, 22574, 33974, 12884, 12944,
                                                 23384, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35108, 0, 3, 22619, 34037, 12944, 13004,
                                                 23474, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35234, 0, 3, 22664, 34100, 13004, 13064,
                                                 23564, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35360, 0, 3, 22709, 34163, 13064, 13124,
                                                 23654, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35486, 0, 3, 22754, 34226, 13124, 13184,
                                                 23744, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35612, 0, 3, 22844, 34289, 13304, 13364,
                                                 23834, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35738, 0, 3, 22889, 34352, 13364, 13424,
                                                 23924, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35864, 0, 3, 22934, 34415, 13424, 13484,
                                                 24014, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35990, 0, 3, 22979, 34478, 13484, 13544,
                                                 24104, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36116, 0, 3, 23024, 34541, 13544, 13604,
                                                 24194, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36242, 0, 3, 23069, 34604, 13604, 13664,
                                                 24284, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36368, 0, 3, 23114, 34667, 13664, 13724,
                                                 24374, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36494, 0, 3, 23294, 34982, 13844, 13944,
                                                 24614, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36704, 0, 3, 23384, 35108, 13944, 14044,
                                                 24764, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36914, 0, 3, 23474, 35234, 14044, 14144,
                                                 24914, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37124, 0, 3, 23564, 35360, 14144, 14244,
                                                 25064, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37334, 0, 3, 23654, 35486, 14244, 14344,
                                                 25214, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37544, 0, 3, 23924, 35864, 14544, 14644,
                                                 25514, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37754, 0, 3, 24014, 35990, 14644, 14744,
                                                 25664, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37964, 0, 3, 24104, 36116, 14744, 14844,
                                                 25814, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38174, 0, 3, 24194, 36242, 14844, 14944,
                                                 25964, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38384, 0, 3, 24284, 36368, 14944, 15044,
                                                 26114, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38594, 0, 3, 34730, 34856, 24464, 36494,
                                                 15244, 15394, 26264, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38909, 0, 3, 34856, 34982, 24614, 36704,
                                                 15394, 15544, 26489, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39224, 0, 3, 34982, 35108, 24764, 36914,
                                                 15544, 15694, 26714, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39539, 0, 3, 35108, 35234, 24914, 37124,
                                                 15694, 15844, 26939, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39854, 0, 3, 35234, 35360, 25064, 37334,
                                                 15844, 15994, 27164, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40169, 0, 3, 35612, 35738, 25364, 37544,
                                                 16294, 16444, 27389, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40484, 0, 3, 35738, 35864, 25514, 37754,
                                                 16444, 16594, 27614, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40799, 0, 3, 35864, 35990, 25664, 37964,
                                                 16594, 16744, 27839, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 41114, 0, 3, 35990, 36116, 25814, 38174,
                                                 16744, 16894, 28064, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 41429, 0, 3, 36116, 36242, 25964, 38384,
                                                 16894, 17044, 28289, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 41744, 0, 3, 36494, 36704, 26489, 39224,
                                                 17344, 17554, 28829, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 42185, 0, 3, 36704, 36914, 26714, 39539,
                                                 17554, 17764, 29144, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 42626, 0, 3, 36914, 37124, 26939, 39854,
                                                 17764, 17974, 29459, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 43067, 0, 3, 37544, 37754, 27614, 40799,
                                                 18394, 18604, 30089, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 43508, 0, 3, 37754, 37964, 27839, 41114,
                                                 18604, 18814, 30404, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 43949, 0, 3, 37964, 38174, 28064, 41429,
                                                 18814, 19024, 30719, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 44390, 0, 3, 38594, 38909, 28514, 41744,
                                                 19444, 19724, 31034, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 44978, 0, 3, 38909, 39224, 28829, 42185,
                                                 19724, 20004, 31454, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 45566, 0, 3, 39224, 39539, 29144, 42626,
                                                 20004, 20284, 31874, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 46154, 0, 3, 40169, 40484, 29774, 43067,
                                                 20844, 21124, 32294, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 46742, 0, 3, 40484, 40799, 30089, 43508,
                                                 21124, 21404, 32714, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 47330, 0, 3, 40799, 41114, 30404, 43949,
                                                 21404, 21684, 33134, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 47918, 3, 22244, 22259, 33575, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 47946, 3, 22259, 22274, 33596, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 47974, 3, 22274, 22289, 33617, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48002, 3, 22289, 22304, 33638, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48030, 3, 22304, 22319, 33659, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48058, 3, 22319, 22334, 33680, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48086, 3, 22364, 22379, 33722, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48114, 3, 22379, 22394, 33743, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48142, 3, 22394, 22409, 33764, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48170, 3, 22409, 22424, 33785, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48198, 3, 22424, 22439, 33806, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 48226, 3, 22439, 22454, 33827, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 48254, 0, 3, 33554, 47918, 33911, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48338, 0, 3, 33575, 47946, 33974, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48422, 0, 3, 33596, 47974, 34037, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48506, 0, 3, 33617, 48002, 34100, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48590, 0, 3, 33638, 48030, 34163, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48674, 0, 3, 33659, 48058, 34226, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48758, 0, 3, 33701, 48086, 34352, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48842, 0, 3, 33722, 48114, 34415, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 48926, 0, 3, 33743, 48142, 34478, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 49010, 0, 3, 33764, 48170, 34541, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 49094, 0, 3, 33785, 48198, 34604, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 49178, 0, 3, 33806, 48226, 34667, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 49262, 0, 3, 33911, 48338, 23204, 23294,
                                                 34982, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 49430, 0, 3, 33974, 48422, 23294, 23384,
                                                 35108, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 49598, 0, 3, 34037, 48506, 23384, 23474,
                                                 35234, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 49766, 0, 3, 34100, 48590, 23474, 23564,
                                                 35360, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 49934, 0, 3, 34163, 48674, 23564, 23654,
                                                 35486, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 50102, 0, 3, 34352, 48842, 23834, 23924,
                                                 35864, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 50270, 0, 3, 34415, 48926, 23924, 24014,
                                                 35990, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 50438, 0, 3, 34478, 49010, 24014, 24104,
                                                 36116, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 50606, 0, 3, 34541, 49094, 24104, 24194,
                                                 36242, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 50774, 0, 3, 34604, 49178, 24194, 24284,
                                                 36368, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 50942, 0, 3, 34982, 49430, 24464, 24614,
                                                 36704, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 51222, 0, 3, 35108, 49598, 24614, 24764,
                                                 36914, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 51502, 0, 3, 35234, 49766, 24764, 24914,
                                                 37124, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 51782, 0, 3, 35360, 49934, 24914, 25064,
                                                 37334, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 52062, 0, 3, 35864, 50270, 25364, 25514,
                                                 37754, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 52342, 0, 3, 35990, 50438, 25514, 25664,
                                                 37964, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 52622, 0, 3, 36116, 50606, 25664, 25814,
                                                 38174, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 52902, 0, 3, 36242, 50774, 25814, 25964,
                                                 38384, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 53182, 0, 3, 49262, 49430, 36704, 51222,
                                                 26264, 26489, 39224, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 53602, 0, 3, 49430, 49598, 36914, 51502,
                                                 26489, 26714, 39539, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 54022, 0, 3, 49598, 49766, 37124, 51782,
                                                 26714, 26939, 39854, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 54442, 0, 3, 50102, 50270, 37754, 52342,
                                                 27389, 27614, 40799, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 54862, 0, 3, 50270, 50438, 37964, 52622,
                                                 27614, 27839, 41114, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 55282, 0, 3, 50438, 50606, 38174, 52902,
                                                 27839, 28064, 41429, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 55702, 0, 3, 50942, 51222, 39224, 53602,
                                                 28514, 28829, 42185, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 56290, 0, 3, 51222, 51502, 39539, 54022,
                                                 28829, 29144, 42626, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 56878, 0, 3, 52062, 52342, 40799, 54862,
                                                 29774, 30089, 43508, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 57466, 0, 3, 52342, 52622, 41114, 55282,
                                                 30089, 30404, 43949, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 58054, 0, 3, 53182, 53602, 42185, 56290,
                                                 31034, 31454, 45566, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 58838, 0, 3, 54442, 54862, 43508, 57466,
                                                 32294, 32714, 47330, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59622, 3, 33554, 33575, 47946, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59658, 3, 33575, 33596, 47974, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59694, 3, 33596, 33617, 48002, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59730, 3, 33617, 33638, 48030, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59766, 3, 33638, 33659, 48058, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59802, 3, 33701, 33722, 48114, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59838, 3, 33722, 33743, 48142, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59874, 3, 33743, 33764, 48170, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59910, 3, 33764, 33785, 48198, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 59946, 3, 33785, 33806, 48226, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 59982, 0, 3, 47918, 59622, 48338, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60090, 0, 3, 47946, 59658, 48422, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60198, 0, 3, 47974, 59694, 48506, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60306, 0, 3, 48002, 59730, 48590, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60414, 0, 3, 48030, 59766, 48674, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60522, 0, 3, 48086, 59802, 48842, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60630, 0, 3, 48114, 59838, 48926, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60738, 0, 3, 48142, 59874, 49010, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60846, 0, 3, 48170, 59910, 49094, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 60954, 0, 3, 48198, 59946, 49178, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 61062, 0, 3, 48254, 59982, 34730, 34856,
                                                 49262, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 61278, 0, 3, 48338, 60090, 34856, 34982,
                                                 49430, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 61494, 0, 3, 48422, 60198, 34982, 35108,
                                                 49598, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 61710, 0, 3, 48506, 60306, 35108, 35234,
                                                 49766, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 61926, 0, 3, 48590, 60414, 35234, 35360,
                                                 49934, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 62142, 0, 3, 48758, 60522, 35612, 35738,
                                                 50102, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 62358, 0, 3, 48842, 60630, 35738, 35864,
                                                 50270, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 62574, 0, 3, 48926, 60738, 35864, 35990,
                                                 50438, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 62790, 0, 3, 49010, 60846, 35990, 36116,
                                                 50606, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 63006, 0, 3, 49094, 60954, 36116, 36242,
                                                 50774, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 63222, 0, 3, 49430, 61494, 36494, 36704,
                                                 51222, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 63582, 0, 3, 49598, 61710, 36704, 36914,
                                                 51502, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 63942, 0, 3, 49766, 61926, 36914, 37124,
                                                 51782, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 64302, 0, 3, 50270, 62574, 37544, 37754,
                                                 52342, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 64662, 0, 3, 50438, 62790, 37754, 37964,
                                                 52622, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 65022, 0, 3, 50606, 63006, 37964, 38174,
                                                 52902, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 65382, 0, 3, 61062, 61278, 50942, 63222,
                                                 38594, 38909, 53182, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 65922, 0, 3, 61278, 61494, 51222, 63582,
                                                 38909, 39224, 53602, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 66462, 0, 3, 61494, 61710, 51502, 63942,
                                                 39224, 39539, 54022, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 67002, 0, 3, 62142, 62358, 52062, 64302,
                                                 40169, 40484, 54442, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 67542, 0, 3, 62358, 62574, 52342, 64662,
                                                 40484, 40799, 54862, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 68082, 0, 3, 62574, 62790, 52622, 65022,
                                                 40799, 41114, 55282, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 68622, 0, 3, 63222, 63582, 53602, 66462,
                                                 41744, 42185, 56290, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 69378, 0, 3, 64302, 64662, 54862, 68082,
                                                 43067, 43508, 57466, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 70134, 0, 3, 65382, 65922, 55702, 68622,
                                                 44390, 44978, 58054, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 71142, 0, 3, 67002, 67542, 56878, 69378,
                                                 46154, 46742, 58838, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 72150, 70134, 2016, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 74166, 73158, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 74166, 15, nmax);

    simdtrf::transform_k_inner(buffer, 74166, 72150, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 195 * nvalues, nvalues, buffer, 74166, 15, nmax);
}

}  // namespace simdt2ceri
