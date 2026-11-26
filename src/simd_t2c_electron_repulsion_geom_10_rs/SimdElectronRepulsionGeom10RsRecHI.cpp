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


#include "SimdElectronRepulsionGeom10RsRecHI.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_hi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_hi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 54409, 50608, 3528, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 12, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 20, 12, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 33, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 7, 8, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 8, 9, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 9, 10, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 10, 11, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 11, 12, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 12, 13, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 13, 14, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 14, 15, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 15, 16, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 16, 17, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 17, 18, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 21, 22, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 22, 23, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 23, 24, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 24, 25, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 25, 26, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 26, 27, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 27, 28, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 28, 29, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 29, 30, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 30, 31, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 31, 32, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 34, 37, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 37, 40, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 40, 43, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 262, 0, 43, 46, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 46, 49, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 282, 0, 49, 52, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 292, 0, 52, 55, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 302, 0, 55, 58, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 58, 61, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 67, 70, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 332, 0, 70, 73, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 342, 0, 73, 76, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 352, 0, 76, 79, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 362, 0, 79, 82, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 372, 0, 82, 85, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 382, 0, 85, 88, 214, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 392, 0, 88, 91, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 402, 0, 91, 94, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 412, 0, 100, 106, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 427, 0, 106, 112, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 442, 0, 112, 118, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 457, 0, 118, 124, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 472, 0, 124, 130, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 487, 0, 130, 136, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 502, 0, 136, 142, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 517, 0, 142, 148, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 532, 0, 148, 154, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 547, 0, 166, 172, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 562, 0, 172, 178, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 577, 0, 178, 184, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 592, 0, 184, 190, 352, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 607, 0, 190, 196, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 622, 0, 196, 202, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 637, 0, 202, 208, 382, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 652, 0, 208, 214, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 667, 0, 214, 220, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 682, 0, 232, 242, 442, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 703, 0, 242, 252, 457, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 724, 0, 252, 262, 472, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 745, 0, 262, 272, 487, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 766, 0, 272, 282, 502, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 787, 0, 282, 292, 517, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 808, 0, 292, 302, 532, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 829, 0, 322, 332, 577, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 850, 0, 332, 342, 592, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 871, 0, 342, 352, 607, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 892, 0, 352, 362, 622, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 913, 0, 362, 372, 637, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 934, 0, 372, 382, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 955, 0, 382, 392, 667, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 976, 0, 412, 427, 682, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1004, 0, 427, 442, 703, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1032, 0, 442, 457, 724, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1060, 0, 457, 472, 745, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1088, 0, 472, 487, 766, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1116, 0, 487, 502, 787, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1144, 0, 502, 517, 808, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1172, 0, 547, 562, 829, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1200, 0, 562, 577, 850, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1228, 0, 577, 592, 871, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1256, 0, 592, 607, 892, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1284, 0, 607, 622, 913, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1312, 0, 622, 637, 934, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1340, 0, 637, 652, 955, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1368, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1371, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1374, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1377, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1380, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1383, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1386, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1389, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1392, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1395, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1398, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1401, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1404, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1407, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1410, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1413, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1416, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1419, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1422, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1425, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1428, 3, 9, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1437, 3, 10, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1446, 3, 11, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1455, 3, 12, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1464, 3, 13, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1473, 3, 14, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1482, 3, 15, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1491, 3, 16, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1500, 3, 17, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1509, 3, 18, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1518, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1527, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1536, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1545, 3, 26, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1554, 3, 27, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1563, 3, 28, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1572, 3, 29, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1581, 3, 30, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1590, 3, 31, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1599, 3, 32, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1608, 0, 3, 37, 1437, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1626, 0, 3, 40, 1446, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1644, 0, 3, 43, 1455, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1662, 0, 3, 46, 1464, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1680, 0, 3, 49, 1473, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1698, 0, 3, 52, 1482, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1716, 0, 3, 55, 1491, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1734, 0, 3, 58, 1500, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1752, 0, 3, 61, 1509, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1770, 0, 3, 70, 1527, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1788, 0, 3, 73, 1536, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1806, 0, 3, 76, 1545, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1824, 0, 3, 79, 1554, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1842, 0, 3, 82, 1563, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1860, 0, 3, 85, 1572, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1878, 0, 3, 88, 1581, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1896, 0, 3, 91, 1590, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1914, 0, 3, 94, 1599, 226, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1932, 0, 3, 112, 1626, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1962, 0, 3, 118, 1644, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1992, 0, 3, 124, 1662, 262, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2022, 0, 3, 130, 1680, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2052, 0, 3, 136, 1698, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2082, 0, 3, 142, 1716, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2112, 0, 3, 148, 1734, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2142, 0, 3, 154, 1752, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2172, 0, 3, 178, 1788, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2202, 0, 3, 184, 1806, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2232, 0, 3, 190, 1824, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2262, 0, 3, 196, 1842, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2292, 0, 3, 202, 1860, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2322, 0, 3, 208, 1878, 382, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2352, 0, 3, 214, 1896, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2382, 0, 3, 220, 1914, 402, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2412, 0, 3, 242, 1962, 442, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2457, 0, 3, 252, 1992, 457, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2502, 0, 3, 262, 2022, 472, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2547, 0, 3, 272, 2052, 487, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2592, 0, 3, 282, 2082, 502, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2637, 0, 3, 292, 2112, 517, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2682, 0, 3, 302, 2142, 532, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2727, 0, 3, 332, 2202, 577, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2772, 0, 3, 342, 2232, 592, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2817, 0, 3, 352, 2262, 607, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2862, 0, 3, 362, 2292, 622, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2907, 0, 3, 372, 2322, 637, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2952, 0, 3, 382, 2352, 652, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2997, 0, 3, 392, 2382, 667, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3042, 0, 3, 442, 2457, 703, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3105, 0, 3, 457, 2502, 724, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3168, 0, 3, 472, 2547, 745, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3231, 0, 3, 487, 2592, 766, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3294, 0, 3, 502, 2637, 787, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3357, 0, 3, 517, 2682, 808, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3420, 0, 3, 577, 2772, 850, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3483, 0, 3, 592, 2817, 871, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3546, 0, 3, 607, 2862, 892, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3609, 0, 3, 622, 2907, 913, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3672, 0, 3, 637, 2952, 934, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3735, 0, 3, 652, 2997, 955, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3798, 0, 3, 703, 3105, 1032, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3882, 0, 3, 724, 3168, 1060, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3966, 0, 3, 745, 3231, 1088, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4050, 0, 3, 766, 3294, 1116, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4134, 0, 3, 787, 3357, 1144, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4218, 0, 3, 850, 3483, 1228, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4302, 0, 3, 871, 3546, 1256, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4386, 0, 3, 892, 3609, 1284, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4470, 0, 3, 913, 3672, 1312, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4554, 0, 3, 934, 3735, 1340, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4638, 3, 9, 10, 1371, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4644, 3, 10, 11, 1374, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4650, 3, 11, 12, 1377, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4656, 3, 12, 13, 1380, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4662, 3, 13, 14, 1383, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4668, 3, 14, 15, 1386, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4674, 3, 15, 16, 1389, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4680, 3, 16, 17, 1392, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4686, 3, 17, 18, 1395, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4692, 3, 23, 24, 1401, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4698, 3, 24, 25, 1404, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4704, 3, 25, 26, 1407, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4710, 3, 26, 27, 1410, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4716, 3, 27, 28, 1413, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4722, 3, 28, 29, 1416, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4728, 3, 29, 30, 1419, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4734, 3, 30, 31, 1422, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4740, 3, 31, 32, 1425, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4746, 0, 3, 1368, 4638, 1437, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4764, 0, 3, 1371, 4644, 1446, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4782, 0, 3, 1374, 4650, 1455, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4800, 0, 3, 1377, 4656, 1464, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4818, 0, 3, 1380, 4662, 1473, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4836, 0, 3, 1383, 4668, 1482, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4854, 0, 3, 1386, 4674, 1491, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4872, 0, 3, 1389, 4680, 1500, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4890, 0, 3, 1392, 4686, 1509, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4908, 0, 3, 1398, 4692, 1527, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4926, 0, 3, 1401, 4698, 1536, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4944, 0, 3, 1404, 4704, 1545, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4962, 0, 3, 1407, 4710, 1554, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4980, 0, 3, 1410, 4716, 1563, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4998, 0, 3, 1413, 4722, 1572, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5016, 0, 3, 1416, 4728, 1581, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5034, 0, 3, 1419, 4734, 1590, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5052, 0, 3, 1422, 4740, 1599, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 5070, 0, 3, 1428, 4746, 100, 106, 1608,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5106, 0, 3, 1437, 4764, 106, 112, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5142, 0, 3, 1446, 4782, 112, 118, 1644,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5178, 0, 3, 1455, 4800, 118, 124, 1662,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5214, 0, 3, 1464, 4818, 124, 130, 1680,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5250, 0, 3, 1473, 4836, 130, 136, 1698,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5286, 0, 3, 1482, 4854, 136, 142, 1716,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5322, 0, 3, 1491, 4872, 142, 148, 1734,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5358, 0, 3, 1500, 4890, 148, 154, 1752,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5394, 0, 3, 1518, 4908, 166, 172, 1770,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5430, 0, 3, 1527, 4926, 172, 178, 1788,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5466, 0, 3, 1536, 4944, 178, 184, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5502, 0, 3, 1545, 4962, 184, 190, 1824,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5538, 0, 3, 1554, 4980, 190, 196, 1842,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5574, 0, 3, 1563, 4998, 196, 202, 1860,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5610, 0, 3, 1572, 5016, 202, 208, 1878,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5646, 0, 3, 1581, 5034, 208, 214, 1896,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5682, 0, 3, 1590, 5052, 214, 220, 1914,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5718, 0, 3, 1626, 5142, 232, 242, 1962,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5778, 0, 3, 1644, 5178, 242, 252, 1992,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5838, 0, 3, 1662, 5214, 252, 262, 2022,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5898, 0, 3, 1680, 5250, 262, 272, 2052,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5958, 0, 3, 1698, 5286, 272, 282, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6018, 0, 3, 1716, 5322, 282, 292, 2112,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6078, 0, 3, 1734, 5358, 292, 302, 2142,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6138, 0, 3, 1788, 5466, 322, 332, 2202,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6198, 0, 3, 1806, 5502, 332, 342, 2232,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6258, 0, 3, 1824, 5538, 342, 352, 2262,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6318, 0, 3, 1842, 5574, 352, 362, 2292,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6378, 0, 3, 1860, 5610, 362, 372, 2322,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6438, 0, 3, 1878, 5646, 372, 382, 2352,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6498, 0, 3, 1896, 5682, 382, 392, 2382,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6558, 0, 3, 5070, 5106, 1932, 5718, 412,
                                                 427, 2412, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6648, 0, 3, 5106, 5142, 1962, 5778, 427,
                                                 442, 2457, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6738, 0, 3, 5142, 5178, 1992, 5838, 442,
                                                 457, 2502, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6828, 0, 3, 5178, 5214, 2022, 5898, 457,
                                                 472, 2547, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6918, 0, 3, 5214, 5250, 2052, 5958, 472,
                                                 487, 2592, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7008, 0, 3, 5250, 5286, 2082, 6018, 487,
                                                 502, 2637, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7098, 0, 3, 5286, 5322, 2112, 6078, 502,
                                                 517, 2682, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7188, 0, 3, 5394, 5430, 2172, 6138, 547,
                                                 562, 2727, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7278, 0, 3, 5430, 5466, 2202, 6198, 562,
                                                 577, 2772, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7368, 0, 3, 5466, 5502, 2232, 6258, 577,
                                                 592, 2817, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7458, 0, 3, 5502, 5538, 2262, 6318, 592,
                                                 607, 2862, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7548, 0, 3, 5538, 5574, 2292, 6378, 607,
                                                 622, 2907, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7638, 0, 3, 5574, 5610, 2322, 6438, 622,
                                                 637, 2952, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7728, 0, 3, 5610, 5646, 2352, 6498, 637,
                                                 652, 2997, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7818, 0, 3, 5718, 5778, 2457, 6738, 682,
                                                 703, 3105, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7944, 0, 3, 5778, 5838, 2502, 6828, 703,
                                                 724, 3168, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8070, 0, 3, 5838, 5898, 2547, 6918, 724,
                                                 745, 3231, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8196, 0, 3, 5898, 5958, 2592, 7008, 745,
                                                 766, 3294, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8322, 0, 3, 5958, 6018, 2637, 7098, 766,
                                                 787, 3357, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8448, 0, 3, 6138, 6198, 2772, 7368, 829,
                                                 850, 3483, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8574, 0, 3, 6198, 6258, 2817, 7458, 850,
                                                 871, 3546, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8700, 0, 3, 6258, 6318, 2862, 7548, 871,
                                                 892, 3609, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8826, 0, 3, 6318, 6378, 2907, 7638, 892,
                                                 913, 3672, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8952, 0, 3, 6378, 6438, 2952, 7728, 913,
                                                 934, 3735, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9078, 0, 3, 6558, 6648, 3042, 7818, 976,
                                                 1004, 3798, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9246, 0, 3, 6648, 6738, 3105, 7944,
                                                 1004, 1032, 3882, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9414, 0, 3, 6738, 6828, 3168, 8070,
                                                 1032, 1060, 3966, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9582, 0, 3, 6828, 6918, 3231, 8196,
                                                 1060, 1088, 4050, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9750, 0, 3, 6918, 7008, 3294, 8322,
                                                 1088, 1116, 4134, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9918, 0, 3, 7188, 7278, 3420, 8448,
                                                 1172, 1200, 4218, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10086, 0, 3, 7278, 7368, 3483, 8574,
                                                 1200, 1228, 4302, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10254, 0, 3, 7368, 7458, 3546, 8700,
                                                 1228, 1256, 4386, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10422, 0, 3, 7458, 7548, 3609, 8826,
                                                 1256, 1284, 4470, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10590, 0, 3, 7548, 7638, 3672, 8952,
                                                 1284, 1312, 4554, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10758, 3, 1368, 1371, 4644, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10768, 3, 1371, 1374, 4650, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10778, 3, 1374, 1377, 4656, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10788, 3, 1377, 1380, 4662, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10798, 3, 1380, 1383, 4668, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10808, 3, 1383, 1386, 4674, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10818, 3, 1386, 1389, 4680, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10828, 3, 1389, 1392, 4686, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10838, 3, 1398, 1401, 4698, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10848, 3, 1401, 1404, 4704, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10858, 3, 1404, 1407, 4710, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10868, 3, 1407, 1410, 4716, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10878, 3, 1410, 1413, 4722, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10888, 3, 1413, 1416, 4728, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10898, 3, 1416, 1419, 4734, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 10908, 3, 1419, 1422, 4740, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 10918, 0, 3, 4638, 10758, 4764, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10948, 0, 3, 4644, 10768, 4782, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 10978, 0, 3, 4650, 10778, 4800, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11008, 0, 3, 4656, 10788, 4818, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11038, 0, 3, 4662, 10798, 4836, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11068, 0, 3, 4668, 10808, 4854, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11098, 0, 3, 4674, 10818, 4872, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11128, 0, 3, 4680, 10828, 4890, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11158, 0, 3, 4692, 10838, 4926, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11188, 0, 3, 4698, 10848, 4944, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11218, 0, 3, 4704, 10858, 4962, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11248, 0, 3, 4710, 10868, 4980, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11278, 0, 3, 4716, 10878, 4998, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11308, 0, 3, 4722, 10888, 5016, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11338, 0, 3, 4728, 10898, 5034, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11368, 0, 3, 4734, 10908, 5052, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 11398, 0, 3, 4764, 10948, 1608, 1626,
                                                 5142, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11458, 0, 3, 4782, 10978, 1626, 1644,
                                                 5178, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11518, 0, 3, 4800, 11008, 1644, 1662,
                                                 5214, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11578, 0, 3, 4818, 11038, 1662, 1680,
                                                 5250, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11638, 0, 3, 4836, 11068, 1680, 1698,
                                                 5286, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11698, 0, 3, 4854, 11098, 1698, 1716,
                                                 5322, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11758, 0, 3, 4872, 11128, 1716, 1734,
                                                 5358, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11818, 0, 3, 4926, 11188, 1770, 1788,
                                                 5466, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11878, 0, 3, 4944, 11218, 1788, 1806,
                                                 5502, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11938, 0, 3, 4962, 11248, 1806, 1824,
                                                 5538, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11998, 0, 3, 4980, 11278, 1824, 1842,
                                                 5574, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12058, 0, 3, 4998, 11308, 1842, 1860,
                                                 5610, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12118, 0, 3, 5016, 11338, 1860, 1878,
                                                 5646, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12178, 0, 3, 5034, 11368, 1878, 1896,
                                                 5682, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12238, 0, 3, 5142, 11458, 1932, 1962,
                                                 5778, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12338, 0, 3, 5178, 11518, 1962, 1992,
                                                 5838, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12438, 0, 3, 5214, 11578, 1992, 2022,
                                                 5898, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12538, 0, 3, 5250, 11638, 2022, 2052,
                                                 5958, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12638, 0, 3, 5286, 11698, 2052, 2082,
                                                 6018, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12738, 0, 3, 5322, 11758, 2082, 2112,
                                                 6078, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12838, 0, 3, 5466, 11878, 2172, 2202,
                                                 6198, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12938, 0, 3, 5502, 11938, 2202, 2232,
                                                 6258, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13038, 0, 3, 5538, 11998, 2232, 2262,
                                                 6318, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13138, 0, 3, 5574, 12058, 2262, 2292,
                                                 6378, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13238, 0, 3, 5610, 12118, 2292, 2322,
                                                 6438, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13338, 0, 3, 5646, 12178, 2322, 2352,
                                                 6498, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13438, 0, 3, 11398, 11458, 5778, 12338,
                                                 2412, 2457, 6738, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13588, 0, 3, 11458, 11518, 5838, 12438,
                                                 2457, 2502, 6828, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13738, 0, 3, 11518, 11578, 5898, 12538,
                                                 2502, 2547, 6918, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13888, 0, 3, 11578, 11638, 5958, 12638,
                                                 2547, 2592, 7008, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14038, 0, 3, 11638, 11698, 6018, 12738,
                                                 2592, 2637, 7098, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14188, 0, 3, 11818, 11878, 6198, 12938,
                                                 2727, 2772, 7368, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14338, 0, 3, 11878, 11938, 6258, 13038,
                                                 2772, 2817, 7458, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14488, 0, 3, 11938, 11998, 6318, 13138,
                                                 2817, 2862, 7548, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14638, 0, 3, 11998, 12058, 6378, 13238,
                                                 2862, 2907, 7638, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14788, 0, 3, 12058, 12118, 6438, 13338,
                                                 2907, 2952, 7728, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14938, 0, 3, 12238, 12338, 6738, 13588,
                                                 3042, 3105, 7944, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15148, 0, 3, 12338, 12438, 6828, 13738,
                                                 3105, 3168, 8070, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15358, 0, 3, 12438, 12538, 6918, 13888,
                                                 3168, 3231, 8196, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15568, 0, 3, 12538, 12638, 7008, 14038,
                                                 3231, 3294, 8322, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15778, 0, 3, 12838, 12938, 7368, 14338,
                                                 3420, 3483, 8574, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15988, 0, 3, 12938, 13038, 7458, 14488,
                                                 3483, 3546, 8700, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16198, 0, 3, 13038, 13138, 7548, 14638,
                                                 3546, 3609, 8826, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16408, 0, 3, 13138, 13238, 7638, 14788,
                                                 3609, 3672, 8952, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16618, 0, 3, 13438, 13588, 7944, 15148,
                                                 3798, 3882, 9414, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16898, 0, 3, 13588, 13738, 8070, 15358,
                                                 3882, 3966, 9582, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17178, 0, 3, 13738, 13888, 8196, 15568,
                                                 3966, 4050, 9750, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17458, 0, 3, 14188, 14338, 8574, 15988,
                                                 4218, 4302, 10254, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17738, 0, 3, 14338, 14488, 8700, 16198,
                                                 4302, 4386, 10422, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18018, 0, 3, 14488, 14638, 8826, 16408,
                                                 4386, 4470, 10590, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18298, 3, 4638, 4644, 10768, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18313, 3, 4644, 4650, 10778, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18328, 3, 4650, 4656, 10788, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18343, 3, 4656, 4662, 10798, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18358, 3, 4662, 4668, 10808, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18373, 3, 4668, 4674, 10818, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18388, 3, 4674, 4680, 10828, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18403, 3, 4692, 4698, 10848, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18418, 3, 4698, 4704, 10858, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18433, 3, 4704, 4710, 10868, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18448, 3, 4710, 4716, 10878, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18463, 3, 4716, 4722, 10888, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18478, 3, 4722, 4728, 10898, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18493, 3, 4728, 4734, 10908, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 18508, 0, 3, 10758, 18298, 10948, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18553, 0, 3, 10768, 18313, 10978, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18598, 0, 3, 10778, 18328, 11008, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18643, 0, 3, 10788, 18343, 11038, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18688, 0, 3, 10798, 18358, 11068, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18733, 0, 3, 10808, 18373, 11098, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18778, 0, 3, 10818, 18388, 11128, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18823, 0, 3, 10838, 18403, 11188, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18868, 0, 3, 10848, 18418, 11218, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18913, 0, 3, 10858, 18433, 11248, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18958, 0, 3, 10868, 18448, 11278, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19003, 0, 3, 10878, 18463, 11308, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19048, 0, 3, 10888, 18478, 11338, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19093, 0, 3, 10898, 18493, 11368, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 19138, 0, 3, 10918, 18508, 5070, 5106,
                                                 11398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19228, 0, 3, 10948, 18553, 5106, 5142,
                                                 11458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19318, 0, 3, 10978, 18598, 5142, 5178,
                                                 11518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19408, 0, 3, 11008, 18643, 5178, 5214,
                                                 11578, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19498, 0, 3, 11038, 18688, 5214, 5250,
                                                 11638, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19588, 0, 3, 11068, 18733, 5250, 5286,
                                                 11698, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19678, 0, 3, 11098, 18778, 5286, 5322,
                                                 11758, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19768, 0, 3, 11158, 18823, 5394, 5430,
                                                 11818, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19858, 0, 3, 11188, 18868, 5430, 5466,
                                                 11878, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19948, 0, 3, 11218, 18913, 5466, 5502,
                                                 11938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20038, 0, 3, 11248, 18958, 5502, 5538,
                                                 11998, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20128, 0, 3, 11278, 19003, 5538, 5574,
                                                 12058, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20218, 0, 3, 11308, 19048, 5574, 5610,
                                                 12118, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20308, 0, 3, 11338, 19093, 5610, 5646,
                                                 12178, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20398, 0, 3, 11458, 19318, 5718, 5778,
                                                 12338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20548, 0, 3, 11518, 19408, 5778, 5838,
                                                 12438, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20698, 0, 3, 11578, 19498, 5838, 5898,
                                                 12538, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20848, 0, 3, 11638, 19588, 5898, 5958,
                                                 12638, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20998, 0, 3, 11698, 19678, 5958, 6018,
                                                 12738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21148, 0, 3, 11878, 19948, 6138, 6198,
                                                 12938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21298, 0, 3, 11938, 20038, 6198, 6258,
                                                 13038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21448, 0, 3, 11998, 20128, 6258, 6318,
                                                 13138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21598, 0, 3, 12058, 20218, 6318, 6378,
                                                 13238, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21748, 0, 3, 12118, 20308, 6378, 6438,
                                                 13338, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21898, 0, 3, 19138, 19228, 12238, 20398,
                                                 6558, 6648, 13438, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22123, 0, 3, 19228, 19318, 12338, 20548,
                                                 6648, 6738, 13588, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22348, 0, 3, 19318, 19408, 12438, 20698,
                                                 6738, 6828, 13738, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22573, 0, 3, 19408, 19498, 12538, 20848,
                                                 6828, 6918, 13888, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22798, 0, 3, 19498, 19588, 12638, 20998,
                                                 6918, 7008, 14038, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23023, 0, 3, 19768, 19858, 12838, 21148,
                                                 7188, 7278, 14188, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23248, 0, 3, 19858, 19948, 12938, 21298,
                                                 7278, 7368, 14338, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23473, 0, 3, 19948, 20038, 13038, 21448,
                                                 7368, 7458, 14488, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23698, 0, 3, 20038, 20128, 13138, 21598,
                                                 7458, 7548, 14638, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23923, 0, 3, 20128, 20218, 13238, 21748,
                                                 7548, 7638, 14788, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24148, 0, 3, 20398, 20548, 13588, 22348,
                                                 7818, 7944, 15148, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24463, 0, 3, 20548, 20698, 13738, 22573,
                                                 7944, 8070, 15358, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24778, 0, 3, 20698, 20848, 13888, 22798,
                                                 8070, 8196, 15568, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 25093, 0, 3, 21148, 21298, 14338, 23473,
                                                 8448, 8574, 15988, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 25408, 0, 3, 21298, 21448, 14488, 23698,
                                                 8574, 8700, 16198, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 25723, 0, 3, 21448, 21598, 14638, 23923,
                                                 8700, 8826, 16408, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26038, 0, 3, 21898, 22123, 14938, 24148,
                                                 9078, 9246, 16618, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26458, 0, 3, 22123, 22348, 15148, 24463,
                                                 9246, 9414, 16898, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26878, 0, 3, 22348, 22573, 15358, 24778,
                                                 9414, 9582, 17178, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 27298, 0, 3, 23023, 23248, 15778, 25093,
                                                 9918, 10086, 17458, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 27718, 0, 3, 23248, 23473, 15988, 25408,
                                                 10086, 10254, 17738, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 28138, 0, 3, 23473, 23698, 16198, 25723,
                                                 10254, 10422, 18018, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28558, 3, 10758, 10768, 18313, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28579, 3, 10768, 10778, 18328, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28600, 3, 10778, 10788, 18343, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28621, 3, 10788, 10798, 18358, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28642, 3, 10798, 10808, 18373, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28663, 3, 10808, 10818, 18388, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28684, 3, 10838, 10848, 18418, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28705, 3, 10848, 10858, 18433, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28726, 3, 10858, 10868, 18448, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28747, 3, 10868, 10878, 18463, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28768, 3, 10878, 10888, 18478, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 28789, 3, 10888, 10898, 18493, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 28810, 0, 3, 18298, 28558, 18553, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 28873, 0, 3, 18313, 28579, 18598, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 28936, 0, 3, 18328, 28600, 18643, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 28999, 0, 3, 18343, 28621, 18688, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29062, 0, 3, 18358, 28642, 18733, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29125, 0, 3, 18373, 28663, 18778, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29188, 0, 3, 18403, 28684, 18868, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29251, 0, 3, 18418, 28705, 18913, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29314, 0, 3, 18433, 28726, 18958, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29377, 0, 3, 18448, 28747, 19003, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29440, 0, 3, 18463, 28768, 19048, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 29503, 0, 3, 18478, 28789, 19093, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 29566, 0, 3, 18553, 28873, 11398, 11458,
                                                 19318, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 29692, 0, 3, 18598, 28936, 11458, 11518,
                                                 19408, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 29818, 0, 3, 18643, 28999, 11518, 11578,
                                                 19498, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 29944, 0, 3, 18688, 29062, 11578, 11638,
                                                 19588, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30070, 0, 3, 18733, 29125, 11638, 11698,
                                                 19678, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30196, 0, 3, 18868, 29251, 11818, 11878,
                                                 19948, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30322, 0, 3, 18913, 29314, 11878, 11938,
                                                 20038, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30448, 0, 3, 18958, 29377, 11938, 11998,
                                                 20128, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30574, 0, 3, 19003, 29440, 11998, 12058,
                                                 20218, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 30700, 0, 3, 19048, 29503, 12058, 12118,
                                                 20308, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 30826, 0, 3, 19318, 29692, 12238, 12338,
                                                 20548, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31036, 0, 3, 19408, 29818, 12338, 12438,
                                                 20698, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31246, 0, 3, 19498, 29944, 12438, 12538,
                                                 20848, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31456, 0, 3, 19588, 30070, 12538, 12638,
                                                 20998, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31666, 0, 3, 19948, 30322, 12838, 12938,
                                                 21298, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 31876, 0, 3, 20038, 30448, 12938, 13038,
                                                 21448, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32086, 0, 3, 20128, 30574, 13038, 13138,
                                                 21598, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32296, 0, 3, 20218, 30700, 13138, 13238,
                                                 21748, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 32506, 0, 3, 29566, 29692, 20548, 31036,
                                                 13438, 13588, 22348, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 32821, 0, 3, 29692, 29818, 20698, 31246,
                                                 13588, 13738, 22573, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33136, 0, 3, 29818, 29944, 20848, 31456,
                                                 13738, 13888, 22798, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33451, 0, 3, 30196, 30322, 21298, 31876,
                                                 14188, 14338, 23473, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33766, 0, 3, 30322, 30448, 21448, 32086,
                                                 14338, 14488, 23698, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34081, 0, 3, 30448, 30574, 21598, 32296,
                                                 14488, 14638, 23923, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34396, 0, 3, 30826, 31036, 22348, 32821,
                                                 14938, 15148, 24463, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34837, 0, 3, 31036, 31246, 22573, 33136,
                                                 15148, 15358, 24778, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35278, 0, 3, 31666, 31876, 23473, 33766,
                                                 15778, 15988, 25408, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35719, 0, 3, 31876, 32086, 23698, 34081,
                                                 15988, 16198, 25723, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 36160, 0, 3, 32506, 32821, 24463, 34837,
                                                 16618, 16898, 26878, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 36748, 0, 3, 33451, 33766, 25408, 35719,
                                                 17458, 17738, 28138, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37336, 3, 18298, 18313, 28579, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37364, 3, 18313, 18328, 28600, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37392, 3, 18328, 18343, 28621, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37420, 3, 18343, 18358, 28642, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37448, 3, 18358, 18373, 28663, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37476, 3, 18403, 18418, 28705, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37504, 3, 18418, 18433, 28726, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37532, 3, 18433, 18448, 28747, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37560, 3, 18448, 18463, 28768, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 37588, 3, 18463, 18478, 28789, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 37616, 0, 3, 28558, 37336, 28873, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 37700, 0, 3, 28579, 37364, 28936, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 37784, 0, 3, 28600, 37392, 28999, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 37868, 0, 3, 28621, 37420, 29062, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 37952, 0, 3, 28642, 37448, 29125, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 38036, 0, 3, 28684, 37476, 29251, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 38120, 0, 3, 28705, 37504, 29314, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 38204, 0, 3, 28726, 37532, 29377, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 38288, 0, 3, 28747, 37560, 29440, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 38372, 0, 3, 28768, 37588, 29503, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 38456, 0, 3, 28810, 37616, 19138, 19228,
                                                 29566, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 38624, 0, 3, 28873, 37700, 19228, 19318,
                                                 29692, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 38792, 0, 3, 28936, 37784, 19318, 19408,
                                                 29818, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 38960, 0, 3, 28999, 37868, 19408, 19498,
                                                 29944, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39128, 0, 3, 29062, 37952, 19498, 19588,
                                                 30070, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39296, 0, 3, 29188, 38036, 19768, 19858,
                                                 30196, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39464, 0, 3, 29251, 38120, 19858, 19948,
                                                 30322, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39632, 0, 3, 29314, 38204, 19948, 20038,
                                                 30448, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39800, 0, 3, 29377, 38288, 20038, 20128,
                                                 30574, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 39968, 0, 3, 29440, 38372, 20128, 20218,
                                                 30700, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40136, 0, 3, 29692, 38792, 20398, 20548,
                                                 31036, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40416, 0, 3, 29818, 38960, 20548, 20698,
                                                 31246, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40696, 0, 3, 29944, 39128, 20698, 20848,
                                                 31456, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40976, 0, 3, 30322, 39632, 21148, 21298,
                                                 31876, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 41256, 0, 3, 30448, 39800, 21298, 21448,
                                                 32086, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 41536, 0, 3, 30574, 39968, 21448, 21598,
                                                 32296, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 41816, 0, 3, 38456, 38624, 30826, 40136,
                                                 21898, 22123, 32506, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 42236, 0, 3, 38624, 38792, 31036, 40416,
                                                 22123, 22348, 32821, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 42656, 0, 3, 38792, 38960, 31246, 40696,
                                                 22348, 22573, 33136, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 43076, 0, 3, 39296, 39464, 31666, 40976,
                                                 23023, 23248, 33451, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 43496, 0, 3, 39464, 39632, 31876, 41256,
                                                 23248, 23473, 33766, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 43916, 0, 3, 39632, 39800, 32086, 41536,
                                                 23473, 23698, 34081, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 44336, 0, 3, 40136, 40416, 32821, 42656,
                                                 24148, 24463, 34837, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 44924, 0, 3, 40976, 41256, 33766, 43916,
                                                 25093, 25408, 35719, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 45512, 0, 3, 41816, 42236, 34396, 44336,
                                                 26038, 26458, 36160, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 46296, 0, 3, 43076, 43496, 35278, 44924,
                                                 27298, 27718, 36748, ncols, alpha, beta, p);

            simdgeo::geom_h_x(buffer, 47080, 43076, 46296, 1, 28, ncols, alpha);

            simdgeo::geom_h_y(buffer, 47668, 43076, 46296, 1, 28, ncols, alpha);

            simdgeo::geom_h_z(buffer, 48256, 43076, 46296, 1, 28, ncols, alpha);

            simdgeo::geom_h_x(buffer, 48844, 41816, 45512, 1, 28, ncols, alpha);

            simdgeo::geom_h_y(buffer, 49432, 41816, 45512, 1, 28, ncols, alpha);

            simdgeo::geom_h_z(buffer, 50020, 41816, 45512, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 50608, 48844, 1764, ncols);

            simdfunc::contract_primitives(buffer, 52372, 47080, 1764, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 54136, 52372, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 54136, 13, nmax);

    simdtrf::transform_i_inner(buffer, 54136, 52960, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 143 * nvalues, nvalues, buffer, 54136, 13, nmax);

    simdtrf::transform_i_inner(buffer, 54136, 53548, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 286 * nvalues, nvalues, buffer, 54136, 13, nmax);

    simdtrf::transform_i_inner(buffer, 54136, 50608, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 429 * nvalues, nvalues, buffer, 54136, 13, nmax);

    simdtrf::transform_i_inner(buffer, 54136, 51196, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 572 * nvalues, nvalues, buffer, 54136, 13, nmax);

    simdtrf::transform_i_inner(buffer, 54136, 51784, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 715 * nvalues, nvalues, buffer, 54136, 13, nmax);
}

}  // namespace simdt2ceri
