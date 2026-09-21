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


#include "SimdElectronRepulsionRsRecHK.hpp"

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
compute_rs_hk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 45305, 43478, 1512, nvalues);

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
                                                9, 10, 11, 12}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 19, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 80, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 83, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 86, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 89, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 92, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 95, 0, 31, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 7, 8, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 8, 9, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 9, 10, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 10, 11, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 11, 12, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 12, 13, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 13, 14, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 14, 15, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 15, 16, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 16, 17, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 20, 21, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 164, 0, 21, 22, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 170, 0, 22, 23, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 176, 0, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 182, 0, 24, 25, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 188, 0, 25, 26, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 194, 0, 26, 27, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 200, 0, 27, 28, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 206, 0, 28, 29, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 212, 0, 29, 30, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 32, 35, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 35, 38, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 38, 41, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 41, 44, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 44, 47, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 47, 50, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 50, 53, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 53, 56, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 56, 59, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 65, 68, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 68, 71, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 71, 74, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 74, 77, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 77, 80, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 80, 83, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 83, 86, 200, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 86, 89, 206, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 89, 92, 212, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 398, 0, 98, 104, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 413, 0, 104, 110, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 428, 0, 110, 116, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 443, 0, 116, 122, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 458, 0, 122, 128, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 473, 0, 128, 134, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 488, 0, 134, 140, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 503, 0, 140, 146, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 518, 0, 158, 164, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 533, 0, 164, 170, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 548, 0, 170, 176, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 563, 0, 176, 182, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 578, 0, 182, 188, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 593, 0, 188, 194, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 608, 0, 194, 200, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 623, 0, 200, 206, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 638, 0, 218, 228, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 659, 0, 228, 238, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 680, 0, 238, 248, 443, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 701, 0, 248, 258, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 722, 0, 258, 268, 473, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 743, 0, 268, 278, 488, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 764, 0, 278, 288, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 785, 0, 308, 318, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 806, 0, 318, 328, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 827, 0, 328, 338, 563, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 848, 0, 338, 348, 578, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 869, 0, 348, 358, 593, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 890, 0, 358, 368, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 911, 0, 368, 378, 623, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 932, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 935, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 938, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 941, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 944, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 947, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 950, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 953, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 956, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 959, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 962, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 965, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 968, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 971, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 974, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 977, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 980, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 983, 3, 31, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 986, 3, 9, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 995, 3, 10, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1004, 3, 11, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1013, 3, 12, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1022, 3, 13, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1031, 3, 14, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1040, 3, 15, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1049, 3, 16, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1058, 3, 17, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1067, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1076, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1085, 3, 24, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1094, 3, 25, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1103, 3, 26, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1112, 3, 27, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1121, 3, 28, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1130, 3, 29, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1139, 3, 30, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1148, 0, 3, 35, 986, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1166, 0, 3, 38, 995, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1184, 0, 3, 41, 1004, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1202, 0, 3, 44, 1013, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1220, 0, 3, 47, 1022, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1238, 0, 3, 50, 1031, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1256, 0, 3, 53, 1040, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1274, 0, 3, 56, 1049, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1292, 0, 3, 59, 1058, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1310, 0, 3, 68, 1067, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1328, 0, 3, 71, 1076, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1346, 0, 3, 74, 1085, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1364, 0, 3, 77, 1094, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1382, 0, 3, 80, 1103, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1400, 0, 3, 83, 1112, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1418, 0, 3, 86, 1121, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1436, 0, 3, 89, 1130, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1454, 0, 3, 92, 1139, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1472, 0, 3, 98, 1148, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1502, 0, 3, 104, 1166, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1532, 0, 3, 110, 1184, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1562, 0, 3, 116, 1202, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1592, 0, 3, 122, 1220, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1622, 0, 3, 128, 1238, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1652, 0, 3, 134, 1256, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1682, 0, 3, 140, 1274, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1712, 0, 3, 146, 1292, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1742, 0, 3, 158, 1310, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1772, 0, 3, 164, 1328, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1802, 0, 3, 170, 1346, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1832, 0, 3, 176, 1364, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1862, 0, 3, 182, 1382, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1892, 0, 3, 188, 1400, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1922, 0, 3, 194, 1418, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1952, 0, 3, 200, 1436, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1982, 0, 3, 206, 1454, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2012, 0, 3, 228, 1532, 413, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2057, 0, 3, 238, 1562, 428, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2102, 0, 3, 248, 1592, 443, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2147, 0, 3, 258, 1622, 458, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2192, 0, 3, 268, 1652, 473, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2237, 0, 3, 278, 1682, 488, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2282, 0, 3, 288, 1712, 503, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2327, 0, 3, 318, 1802, 533, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2372, 0, 3, 328, 1832, 548, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2417, 0, 3, 338, 1862, 563, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2462, 0, 3, 348, 1892, 578, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2507, 0, 3, 358, 1922, 593, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2552, 0, 3, 368, 1952, 608, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2597, 0, 3, 378, 1982, 623, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2642, 0, 3, 398, 2012, 638, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2705, 0, 3, 413, 2057, 659, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2768, 0, 3, 428, 2102, 680, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2831, 0, 3, 443, 2147, 701, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2894, 0, 3, 458, 2192, 722, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2957, 0, 3, 473, 2237, 743, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3020, 0, 3, 488, 2282, 764, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3083, 0, 3, 518, 2327, 785, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3146, 0, 3, 533, 2372, 806, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3209, 0, 3, 548, 2417, 827, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3272, 0, 3, 563, 2462, 848, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3335, 0, 3, 578, 2507, 869, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3398, 0, 3, 593, 2552, 890, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3461, 0, 3, 608, 2597, 911, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3524, 3, 9, 10, 935, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 3530, 3, 10, 11, 938, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3536, 3, 11, 12, 941, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3542, 3, 12, 13, 944, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3548, 3, 13, 14, 947, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3554, 3, 14, 15, 950, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3560, 3, 15, 16, 953, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3566, 3, 16, 17, 956, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3572, 3, 22, 23, 962, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3578, 3, 23, 24, 965, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3584, 3, 24, 25, 968, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3590, 3, 25, 26, 971, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3596, 3, 26, 27, 974, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3602, 3, 27, 28, 977, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3608, 3, 28, 29, 980, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3614, 3, 29, 30, 983, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3620, 0, 3, 932, 3524, 995, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3638, 0, 3, 935, 3530, 1004, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3656, 0, 3, 938, 3536, 1013, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3674, 0, 3, 941, 3542, 1022, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3692, 0, 3, 944, 3548, 1031, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3710, 0, 3, 947, 3554, 1040, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3728, 0, 3, 950, 3560, 1049, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3746, 0, 3, 953, 3566, 1058, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3764, 0, 3, 959, 3572, 1076, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3782, 0, 3, 962, 3578, 1085, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3800, 0, 3, 965, 3584, 1094, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3818, 0, 3, 968, 3590, 1103, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3836, 0, 3, 971, 3596, 1112, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3854, 0, 3, 974, 3602, 1121, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3872, 0, 3, 977, 3608, 1130, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3890, 0, 3, 980, 3614, 1139, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3908, 0, 3, 986, 3620, 98, 104, 1166,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3944, 0, 3, 995, 3638, 104, 110, 1184,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3980, 0, 3, 1004, 3656, 110, 116, 1202,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4016, 0, 3, 1013, 3674, 116, 122, 1220,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4052, 0, 3, 1022, 3692, 122, 128, 1238,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4088, 0, 3, 1031, 3710, 128, 134, 1256,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4124, 0, 3, 1040, 3728, 134, 140, 1274,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4160, 0, 3, 1049, 3746, 140, 146, 1292,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4196, 0, 3, 1067, 3764, 158, 164, 1328,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4232, 0, 3, 1076, 3782, 164, 170, 1346,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4268, 0, 3, 1085, 3800, 170, 176, 1364,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4304, 0, 3, 1094, 3818, 176, 182, 1382,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4340, 0, 3, 1103, 3836, 182, 188, 1400,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4376, 0, 3, 1112, 3854, 188, 194, 1418,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4412, 0, 3, 1121, 3872, 194, 200, 1436,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4448, 0, 3, 1130, 3890, 200, 206, 1454,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4484, 0, 3, 1166, 3944, 218, 228, 1532,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4544, 0, 3, 1184, 3980, 228, 238, 1562,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4604, 0, 3, 1202, 4016, 238, 248, 1592,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4664, 0, 3, 1220, 4052, 248, 258, 1622,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4724, 0, 3, 1238, 4088, 258, 268, 1652,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4784, 0, 3, 1256, 4124, 268, 278, 1682,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4844, 0, 3, 1274, 4160, 278, 288, 1712,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4904, 0, 3, 1328, 4232, 308, 318, 1802,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4964, 0, 3, 1346, 4268, 318, 328, 1832,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5024, 0, 3, 1364, 4304, 328, 338, 1862,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5084, 0, 3, 1382, 4340, 338, 348, 1892,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5144, 0, 3, 1400, 4376, 348, 358, 1922,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5204, 0, 3, 1418, 4412, 358, 368, 1952,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5264, 0, 3, 1436, 4448, 368, 378, 1982,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5324, 0, 3, 3908, 3944, 1532, 4544, 398,
                                                 413, 2057, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5414, 0, 3, 3944, 3980, 1562, 4604, 413,
                                                 428, 2102, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5504, 0, 3, 3980, 4016, 1592, 4664, 428,
                                                 443, 2147, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5594, 0, 3, 4016, 4052, 1622, 4724, 443,
                                                 458, 2192, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5684, 0, 3, 4052, 4088, 1652, 4784, 458,
                                                 473, 2237, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5774, 0, 3, 4088, 4124, 1682, 4844, 473,
                                                 488, 2282, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5864, 0, 3, 4196, 4232, 1802, 4964, 518,
                                                 533, 2372, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5954, 0, 3, 4232, 4268, 1832, 5024, 533,
                                                 548, 2417, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6044, 0, 3, 4268, 4304, 1862, 5084, 548,
                                                 563, 2462, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6134, 0, 3, 4304, 4340, 1892, 5144, 563,
                                                 578, 2507, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6224, 0, 3, 4340, 4376, 1922, 5204, 578,
                                                 593, 2552, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6314, 0, 3, 4376, 4412, 1952, 5264, 593,
                                                 608, 2597, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6404, 0, 3, 4484, 4544, 2057, 5414, 638,
                                                 659, 2768, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6530, 0, 3, 4544, 4604, 2102, 5504, 659,
                                                 680, 2831, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6656, 0, 3, 4604, 4664, 2147, 5594, 680,
                                                 701, 2894, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6782, 0, 3, 4664, 4724, 2192, 5684, 701,
                                                 722, 2957, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6908, 0, 3, 4724, 4784, 2237, 5774, 722,
                                                 743, 3020, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7034, 0, 3, 4904, 4964, 2372, 5954, 785,
                                                 806, 3209, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7160, 0, 3, 4964, 5024, 2417, 6044, 806,
                                                 827, 3272, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7286, 0, 3, 5024, 5084, 2462, 6134, 827,
                                                 848, 3335, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7412, 0, 3, 5084, 5144, 2507, 6224, 848,
                                                 869, 3398, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7538, 0, 3, 5144, 5204, 2552, 6314, 869,
                                                 890, 3461, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7664, 3, 932, 935, 3530, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7674, 3, 935, 938, 3536, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7684, 3, 938, 941, 3542, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7694, 3, 941, 944, 3548, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7704, 3, 944, 947, 3554, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7714, 3, 947, 950, 3560, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7724, 3, 950, 953, 3566, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7734, 3, 959, 962, 3578, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7744, 3, 962, 965, 3584, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7754, 3, 965, 968, 3590, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7764, 3, 968, 971, 3596, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7774, 3, 971, 974, 3602, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7784, 3, 974, 977, 3608, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7794, 3, 977, 980, 3614, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 7804, 0, 3, 3524, 7664, 3638, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7834, 0, 3, 3530, 7674, 3656, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7864, 0, 3, 3536, 7684, 3674, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7894, 0, 3, 3542, 7694, 3692, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7924, 0, 3, 3548, 7704, 3710, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7954, 0, 3, 3554, 7714, 3728, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7984, 0, 3, 3560, 7724, 3746, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8014, 0, 3, 3572, 7734, 3782, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8044, 0, 3, 3578, 7744, 3800, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8074, 0, 3, 3584, 7754, 3818, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8104, 0, 3, 3590, 7764, 3836, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8134, 0, 3, 3596, 7774, 3854, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8164, 0, 3, 3602, 7784, 3872, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8194, 0, 3, 3608, 7794, 3890, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 8224, 0, 3, 3620, 7804, 1148, 1166,
                                                 3944, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8284, 0, 3, 3638, 7834, 1166, 1184,
                                                 3980, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8344, 0, 3, 3656, 7864, 1184, 1202,
                                                 4016, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8404, 0, 3, 3674, 7894, 1202, 1220,
                                                 4052, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8464, 0, 3, 3692, 7924, 1220, 1238,
                                                 4088, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8524, 0, 3, 3710, 7954, 1238, 1256,
                                                 4124, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8584, 0, 3, 3728, 7984, 1256, 1274,
                                                 4160, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8644, 0, 3, 3764, 8014, 1310, 1328,
                                                 4232, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8704, 0, 3, 3782, 8044, 1328, 1346,
                                                 4268, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8764, 0, 3, 3800, 8074, 1346, 1364,
                                                 4304, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8824, 0, 3, 3818, 8104, 1364, 1382,
                                                 4340, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8884, 0, 3, 3836, 8134, 1382, 1400,
                                                 4376, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8944, 0, 3, 3854, 8164, 1400, 1418,
                                                 4412, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9004, 0, 3, 3872, 8194, 1418, 1436,
                                                 4448, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9064, 0, 3, 3908, 8224, 1472, 1502,
                                                 4484, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9164, 0, 3, 3944, 8284, 1502, 1532,
                                                 4544, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9264, 0, 3, 3980, 8344, 1532, 1562,
                                                 4604, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9364, 0, 3, 4016, 8404, 1562, 1592,
                                                 4664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9464, 0, 3, 4052, 8464, 1592, 1622,
                                                 4724, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9564, 0, 3, 4088, 8524, 1622, 1652,
                                                 4784, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9664, 0, 3, 4124, 8584, 1652, 1682,
                                                 4844, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9764, 0, 3, 4196, 8644, 1742, 1772,
                                                 4904, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9864, 0, 3, 4232, 8704, 1772, 1802,
                                                 4964, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9964, 0, 3, 4268, 8764, 1802, 1832,
                                                 5024, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10064, 0, 3, 4304, 8824, 1832, 1862,
                                                 5084, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10164, 0, 3, 4340, 8884, 1862, 1892,
                                                 5144, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10264, 0, 3, 4376, 8944, 1892, 1922,
                                                 5204, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10364, 0, 3, 4412, 9004, 1922, 1952,
                                                 5264, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10464, 0, 3, 8224, 8284, 4544, 9264,
                                                 2012, 2057, 5414, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10614, 0, 3, 8284, 8344, 4604, 9364,
                                                 2057, 2102, 5504, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10764, 0, 3, 8344, 8404, 4664, 9464,
                                                 2102, 2147, 5594, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10914, 0, 3, 8404, 8464, 4724, 9564,
                                                 2147, 2192, 5684, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11064, 0, 3, 8464, 8524, 4784, 9664,
                                                 2192, 2237, 5774, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11214, 0, 3, 8644, 8704, 4964, 9964,
                                                 2327, 2372, 5954, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11364, 0, 3, 8704, 8764, 5024, 10064,
                                                 2372, 2417, 6044, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11514, 0, 3, 8764, 8824, 5084, 10164,
                                                 2417, 2462, 6134, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11664, 0, 3, 8824, 8884, 5144, 10264,
                                                 2462, 2507, 6224, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 11814, 0, 3, 8884, 8944, 5204, 10364,
                                                 2507, 2552, 6314, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11964, 0, 3, 9064, 9164, 5324, 10464,
                                                 2642, 2705, 6404, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12174, 0, 3, 9164, 9264, 5414, 10614,
                                                 2705, 2768, 6530, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12384, 0, 3, 9264, 9364, 5504, 10764,
                                                 2768, 2831, 6656, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12594, 0, 3, 9364, 9464, 5594, 10914,
                                                 2831, 2894, 6782, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 12804, 0, 3, 9464, 9564, 5684, 11064,
                                                 2894, 2957, 6908, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13014, 0, 3, 9764, 9864, 5864, 11214,
                                                 3083, 3146, 7034, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13224, 0, 3, 9864, 9964, 5954, 11364,
                                                 3146, 3209, 7160, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13434, 0, 3, 9964, 10064, 6044, 11514,
                                                 3209, 3272, 7286, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13644, 0, 3, 10064, 10164, 6134, 11664,
                                                 3272, 3335, 7412, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13854, 0, 3, 10164, 10264, 6224, 11814,
                                                 3335, 3398, 7538, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14064, 3, 3524, 3530, 7674, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14079, 3, 3530, 3536, 7684, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14094, 3, 3536, 3542, 7694, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14109, 3, 3542, 3548, 7704, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14124, 3, 3548, 3554, 7714, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14139, 3, 3554, 3560, 7724, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14154, 3, 3572, 3578, 7744, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14169, 3, 3578, 3584, 7754, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14184, 3, 3584, 3590, 7764, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14199, 3, 3590, 3596, 7774, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14214, 3, 3596, 3602, 7784, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 14229, 3, 3602, 3608, 7794, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 14244, 0, 3, 7664, 14064, 7834, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14289, 0, 3, 7674, 14079, 7864, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14334, 0, 3, 7684, 14094, 7894, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14379, 0, 3, 7694, 14109, 7924, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14424, 0, 3, 7704, 14124, 7954, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14469, 0, 3, 7714, 14139, 7984, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14514, 0, 3, 7734, 14154, 8044, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14559, 0, 3, 7744, 14169, 8074, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14604, 0, 3, 7754, 14184, 8104, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14649, 0, 3, 7764, 14199, 8134, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14694, 0, 3, 7774, 14214, 8164, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 14739, 0, 3, 7784, 14229, 8194, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 14784, 0, 3, 7804, 14244, 3908, 3944,
                                                 8284, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14874, 0, 3, 7834, 14289, 3944, 3980,
                                                 8344, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14964, 0, 3, 7864, 14334, 3980, 4016,
                                                 8404, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15054, 0, 3, 7894, 14379, 4016, 4052,
                                                 8464, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15144, 0, 3, 7924, 14424, 4052, 4088,
                                                 8524, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15234, 0, 3, 7954, 14469, 4088, 4124,
                                                 8584, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15324, 0, 3, 8014, 14514, 4196, 4232,
                                                 8704, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15414, 0, 3, 8044, 14559, 4232, 4268,
                                                 8764, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15504, 0, 3, 8074, 14604, 4268, 4304,
                                                 8824, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15594, 0, 3, 8104, 14649, 4304, 4340,
                                                 8884, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15684, 0, 3, 8134, 14694, 4340, 4376,
                                                 8944, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 15774, 0, 3, 8164, 14739, 4376, 4412,
                                                 9004, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 15864, 0, 3, 8284, 14874, 4484, 4544,
                                                 9264, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16014, 0, 3, 8344, 14964, 4544, 4604,
                                                 9364, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16164, 0, 3, 8404, 15054, 4604, 4664,
                                                 9464, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16314, 0, 3, 8464, 15144, 4664, 4724,
                                                 9564, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16464, 0, 3, 8524, 15234, 4724, 4784,
                                                 9664, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16614, 0, 3, 8704, 15414, 4904, 4964,
                                                 9964, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16764, 0, 3, 8764, 15504, 4964, 5024,
                                                 10064, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 16914, 0, 3, 8824, 15594, 5024, 5084,
                                                 10164, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17064, 0, 3, 8884, 15684, 5084, 5144,
                                                 10264, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 17214, 0, 3, 8944, 15774, 5144, 5204,
                                                 10364, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 17364, 0, 3, 14784, 14874, 9264, 16014,
                                                 5324, 5414, 10614, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 17589, 0, 3, 14874, 14964, 9364, 16164,
                                                 5414, 5504, 10764, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 17814, 0, 3, 14964, 15054, 9464, 16314,
                                                 5504, 5594, 10914, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18039, 0, 3, 15054, 15144, 9564, 16464,
                                                 5594, 5684, 11064, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18264, 0, 3, 15324, 15414, 9964, 16764,
                                                 5864, 5954, 11364, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18489, 0, 3, 15414, 15504, 10064, 16914,
                                                 5954, 6044, 11514, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18714, 0, 3, 15504, 15594, 10164, 17064,
                                                 6044, 6134, 11664, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 18939, 0, 3, 15594, 15684, 10264, 17214,
                                                 6134, 6224, 11814, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 19164, 0, 3, 15864, 16014, 10614, 17589,
                                                 6404, 6530, 12384, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 19479, 0, 3, 16014, 16164, 10764, 17814,
                                                 6530, 6656, 12594, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 19794, 0, 3, 16164, 16314, 10914, 18039,
                                                 6656, 6782, 12804, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20109, 0, 3, 16614, 16764, 11364, 18489,
                                                 7034, 7160, 13434, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20424, 0, 3, 16764, 16914, 11514, 18714,
                                                 7160, 7286, 13644, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 20739, 0, 3, 16914, 17064, 11664, 18939,
                                                 7286, 7412, 13854, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21054, 3, 7664, 7674, 14079, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21075, 3, 7674, 7684, 14094, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21096, 3, 7684, 7694, 14109, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21117, 3, 7694, 7704, 14124, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21138, 3, 7704, 7714, 14139, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21159, 3, 7734, 7744, 14169, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21180, 3, 7744, 7754, 14184, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21201, 3, 7754, 7764, 14199, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21222, 3, 7764, 7774, 14214, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21243, 3, 7774, 7784, 14229, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 21264, 0, 3, 14064, 21054, 14289, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21327, 0, 3, 14079, 21075, 14334, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21390, 0, 3, 14094, 21096, 14379, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21453, 0, 3, 14109, 21117, 14424, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21516, 0, 3, 14124, 21138, 14469, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21579, 0, 3, 14154, 21159, 14559, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21642, 0, 3, 14169, 21180, 14604, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21705, 0, 3, 14184, 21201, 14649, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21768, 0, 3, 14199, 21222, 14694, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21831, 0, 3, 14214, 21243, 14739, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 21894, 0, 3, 14244, 21264, 8224, 8284,
                                                 14874, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22020, 0, 3, 14289, 21327, 8284, 8344,
                                                 14964, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22146, 0, 3, 14334, 21390, 8344, 8404,
                                                 15054, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22272, 0, 3, 14379, 21453, 8404, 8464,
                                                 15144, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22398, 0, 3, 14424, 21516, 8464, 8524,
                                                 15234, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22524, 0, 3, 14514, 21579, 8644, 8704,
                                                 15414, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22650, 0, 3, 14559, 21642, 8704, 8764,
                                                 15504, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22776, 0, 3, 14604, 21705, 8764, 8824,
                                                 15594, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22902, 0, 3, 14649, 21768, 8824, 8884,
                                                 15684, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 23028, 0, 3, 14694, 21831, 8884, 8944,
                                                 15774, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23154, 0, 3, 14784, 21894, 9064, 9164,
                                                 15864, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23364, 0, 3, 14874, 22020, 9164, 9264,
                                                 16014, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23574, 0, 3, 14964, 22146, 9264, 9364,
                                                 16164, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23784, 0, 3, 15054, 22272, 9364, 9464,
                                                 16314, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23994, 0, 3, 15144, 22398, 9464, 9564,
                                                 16464, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24204, 0, 3, 15324, 22524, 9764, 9864,
                                                 16614, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24414, 0, 3, 15414, 22650, 9864, 9964,
                                                 16764, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24624, 0, 3, 15504, 22776, 9964, 10064,
                                                 16914, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 24834, 0, 3, 15594, 22902, 10064, 10164,
                                                 17064, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 25044, 0, 3, 15684, 23028, 10164, 10264,
                                                 17214, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25254, 0, 3, 21894, 22020, 16014, 23574,
                                                 10464, 10614, 17589, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25569, 0, 3, 22020, 22146, 16164, 23784,
                                                 10614, 10764, 17814, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 25884, 0, 3, 22146, 22272, 16314, 23994,
                                                 10764, 10914, 18039, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26199, 0, 3, 22524, 22650, 16764, 24624,
                                                 11214, 11364, 18489, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26514, 0, 3, 22650, 22776, 16914, 24834,
                                                 11364, 11514, 18714, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 26829, 0, 3, 22776, 22902, 17064, 25044,
                                                 11514, 11664, 18939, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 27144, 0, 3, 23154, 23364, 17364, 25254,
                                                 11964, 12174, 19164, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 27585, 0, 3, 23364, 23574, 17589, 25569,
                                                 12174, 12384, 19479, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 28026, 0, 3, 23574, 23784, 17814, 25884,
                                                 12384, 12594, 19794, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 28467, 0, 3, 24204, 24414, 18264, 26199,
                                                 13014, 13224, 20109, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 28908, 0, 3, 24414, 24624, 18489, 26514,
                                                 13224, 13434, 20424, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 29349, 0, 3, 24624, 24834, 18714, 26829,
                                                 13434, 13644, 20739, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29790, 3, 14064, 14079, 21075, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29818, 3, 14079, 14094, 21096, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29846, 3, 14094, 14109, 21117, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29874, 3, 14109, 14124, 21138, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29902, 3, 14154, 14169, 21180, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29930, 3, 14169, 14184, 21201, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29958, 3, 14184, 14199, 21222, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 29986, 3, 14199, 14214, 21243, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 30014, 0, 3, 21054, 29790, 21327, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30098, 0, 3, 21075, 29818, 21390, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30182, 0, 3, 21096, 29846, 21453, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30266, 0, 3, 21117, 29874, 21516, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30350, 0, 3, 21159, 29902, 21642, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30434, 0, 3, 21180, 29930, 21705, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30518, 0, 3, 21201, 29958, 21768, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 30602, 0, 3, 21222, 29986, 21831, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 30686, 0, 3, 21264, 30014, 14784, 14874,
                                                 22020, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 30854, 0, 3, 21327, 30098, 14874, 14964,
                                                 22146, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31022, 0, 3, 21390, 30182, 14964, 15054,
                                                 22272, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31190, 0, 3, 21453, 30266, 15054, 15144,
                                                 22398, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31358, 0, 3, 21579, 30350, 15324, 15414,
                                                 22650, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31526, 0, 3, 21642, 30434, 15414, 15504,
                                                 22776, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31694, 0, 3, 21705, 30518, 15504, 15594,
                                                 22902, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 31862, 0, 3, 21768, 30602, 15594, 15684,
                                                 23028, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 32030, 0, 3, 22020, 30854, 15864, 16014,
                                                 23574, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 32310, 0, 3, 22146, 31022, 16014, 16164,
                                                 23784, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 32590, 0, 3, 22272, 31190, 16164, 16314,
                                                 23994, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 32870, 0, 3, 22650, 31526, 16614, 16764,
                                                 24624, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 33150, 0, 3, 22776, 31694, 16764, 16914,
                                                 24834, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 33430, 0, 3, 22902, 31862, 16914, 17064,
                                                 25044, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 33710, 0, 3, 30686, 30854, 23574, 32310,
                                                 17364, 17589, 25569, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 34130, 0, 3, 30854, 31022, 23784, 32590,
                                                 17589, 17814, 25884, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 34550, 0, 3, 31358, 31526, 24624, 33150,
                                                 18264, 18489, 26514, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 34970, 0, 3, 31526, 31694, 24834, 33430,
                                                 18489, 18714, 26829, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 35390, 0, 3, 32030, 32310, 25569, 34130,
                                                 19164, 19479, 28026, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 35978, 0, 3, 32870, 33150, 26514, 34970,
                                                 20109, 20424, 29349, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36566, 3, 21054, 21075, 29818, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36602, 3, 21075, 21096, 29846, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36638, 3, 21096, 21117, 29874, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36674, 3, 21159, 21180, 29930, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36710, 3, 21180, 21201, 29958, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 36746, 3, 21201, 21222, 29986, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 36782, 0, 3, 29790, 36566, 30098, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 36890, 0, 3, 29818, 36602, 30182, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 36998, 0, 3, 29846, 36638, 30266, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 37106, 0, 3, 29902, 36674, 30434, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 37214, 0, 3, 29930, 36710, 30518, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 37322, 0, 3, 29958, 36746, 30602, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 37430, 0, 3, 30014, 36782, 21894, 22020,
                                                 30854, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 37646, 0, 3, 30098, 36890, 22020, 22146,
                                                 31022, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 37862, 0, 3, 30182, 36998, 22146, 22272,
                                                 31190, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 38078, 0, 3, 30350, 37106, 22524, 22650,
                                                 31526, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 38294, 0, 3, 30434, 37214, 22650, 22776,
                                                 31694, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 38510, 0, 3, 30518, 37322, 22776, 22902,
                                                 31862, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 38726, 0, 3, 30686, 37430, 23154, 23364,
                                                 32030, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 39086, 0, 3, 30854, 37646, 23364, 23574,
                                                 32310, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 39446, 0, 3, 31022, 37862, 23574, 23784,
                                                 32590, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 39806, 0, 3, 31358, 38078, 24204, 24414,
                                                 32870, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 40166, 0, 3, 31526, 38294, 24414, 24624,
                                                 33150, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 40526, 0, 3, 31694, 38510, 24624, 24834,
                                                 33430, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 40886, 0, 3, 37430, 37646, 32310, 39446,
                                                 25254, 25569, 34130, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 41426, 0, 3, 38078, 38294, 33150, 40526,
                                                 26199, 26514, 34970, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 41966, 0, 3, 38726, 39086, 33710, 40886,
                                                 27144, 27585, 35390, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 42722, 0, 3, 39806, 40166, 34550, 41426,
                                                 28467, 28908, 35978, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 43478, 41966, 1512, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 44990, 44234, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 44990, 15, nmax);

    simdtrf::transform_k_inner(buffer, 44990, 43478, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 165 * nvalues, nvalues, buffer, 44990, 15, nmax);
}

}  // namespace simdt2ceri
