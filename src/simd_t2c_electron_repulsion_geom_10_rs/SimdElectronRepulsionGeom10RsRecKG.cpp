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


#include "SimdElectronRepulsionGeom10RsRecKG.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_kg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_kg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 45532, 41968, 3240, nvalues);

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

            compute_prim_ks_electron_repulsion_0(buffer, 1368, 0, 682, 703, 1032, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1404, 0, 703, 724, 1060, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1440, 0, 724, 745, 1088, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1476, 0, 745, 766, 1116, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1512, 0, 766, 787, 1144, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1548, 0, 829, 850, 1228, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1584, 0, 850, 871, 1256, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1620, 0, 871, 892, 1284, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1656, 0, 892, 913, 1312, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1692, 0, 913, 934, 1340, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1728, 0, 976, 1004, 1368, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1773, 0, 1004, 1032, 1404, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1818, 0, 1032, 1060, 1440, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1863, 0, 1060, 1088, 1476, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1908, 0, 1088, 1116, 1512, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1953, 0, 1172, 1200, 1548, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1998, 0, 1200, 1228, 1584, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2043, 0, 1228, 1256, 1620, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2088, 0, 1256, 1284, 1656, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2133, 0, 1284, 1312, 1692, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2178, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2181, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2184, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2187, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2190, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2193, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2196, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2199, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2202, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2205, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2208, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2211, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2214, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2217, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2220, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2223, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2226, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2229, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2232, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2235, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2238, 3, 9, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2247, 3, 10, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2256, 3, 11, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2265, 3, 12, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2274, 3, 13, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2283, 3, 14, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2292, 3, 15, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2301, 3, 16, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2310, 3, 17, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2319, 3, 18, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2328, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2337, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2346, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2355, 3, 26, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2364, 3, 27, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2373, 3, 28, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2382, 3, 29, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2391, 3, 30, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2400, 3, 31, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2409, 3, 32, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2418, 0, 3, 37, 2247, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2436, 0, 3, 40, 2256, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2454, 0, 3, 43, 2265, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2472, 0, 3, 46, 2274, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2490, 0, 3, 49, 2283, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2508, 0, 3, 52, 2292, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2526, 0, 3, 55, 2301, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2544, 0, 3, 58, 2310, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2562, 0, 3, 61, 2319, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2580, 0, 3, 70, 2337, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2598, 0, 3, 73, 2346, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2616, 0, 3, 76, 2355, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2634, 0, 3, 79, 2364, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2652, 0, 3, 82, 2373, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2670, 0, 3, 85, 2382, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2688, 0, 3, 88, 2391, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2706, 0, 3, 91, 2400, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2724, 0, 3, 94, 2409, 226, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2742, 0, 3, 112, 2436, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2772, 0, 3, 118, 2454, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2802, 0, 3, 124, 2472, 262, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2832, 0, 3, 130, 2490, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2862, 0, 3, 136, 2508, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2892, 0, 3, 142, 2526, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2922, 0, 3, 148, 2544, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2952, 0, 3, 154, 2562, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2982, 0, 3, 178, 2598, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3012, 0, 3, 184, 2616, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3042, 0, 3, 190, 2634, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3072, 0, 3, 196, 2652, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3102, 0, 3, 202, 2670, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3132, 0, 3, 208, 2688, 382, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3162, 0, 3, 214, 2706, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3192, 0, 3, 220, 2724, 402, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3222, 0, 3, 242, 2772, 442, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3267, 0, 3, 252, 2802, 457, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3312, 0, 3, 262, 2832, 472, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3357, 0, 3, 272, 2862, 487, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3402, 0, 3, 282, 2892, 502, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3447, 0, 3, 292, 2922, 517, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3492, 0, 3, 302, 2952, 532, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3537, 0, 3, 332, 3012, 577, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3582, 0, 3, 342, 3042, 592, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3627, 0, 3, 352, 3072, 607, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3672, 0, 3, 362, 3102, 622, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3717, 0, 3, 372, 3132, 637, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3762, 0, 3, 382, 3162, 652, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3807, 0, 3, 392, 3192, 667, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3852, 0, 3, 442, 3267, 703, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3915, 0, 3, 457, 3312, 724, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3978, 0, 3, 472, 3357, 745, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4041, 0, 3, 487, 3402, 766, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4104, 0, 3, 502, 3447, 787, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4167, 0, 3, 517, 3492, 808, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4230, 0, 3, 577, 3582, 850, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4293, 0, 3, 592, 3627, 871, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4356, 0, 3, 607, 3672, 892, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4419, 0, 3, 622, 3717, 913, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4482, 0, 3, 637, 3762, 934, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4545, 0, 3, 652, 3807, 955, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4608, 0, 3, 703, 3915, 1032, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4692, 0, 3, 724, 3978, 1060, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4776, 0, 3, 745, 4041, 1088, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4860, 0, 3, 766, 4104, 1116, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4944, 0, 3, 787, 4167, 1144, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5028, 0, 3, 850, 4293, 1228, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5112, 0, 3, 871, 4356, 1256, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5196, 0, 3, 892, 4419, 1284, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5280, 0, 3, 913, 4482, 1312, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5364, 0, 3, 934, 4545, 1340, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5448, 0, 3, 1032, 4692, 1404, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5556, 0, 3, 1060, 4776, 1440, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5664, 0, 3, 1088, 4860, 1476, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5772, 0, 3, 1116, 4944, 1512, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5880, 0, 3, 1228, 5112, 1584, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5988, 0, 3, 1256, 5196, 1620, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6096, 0, 3, 1284, 5280, 1656, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6204, 0, 3, 1312, 5364, 1692, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6312, 0, 3, 1404, 5556, 1818, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6447, 0, 3, 1440, 5664, 1863, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6582, 0, 3, 1476, 5772, 1908, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6717, 0, 3, 1584, 5988, 2043, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6852, 0, 3, 1620, 6096, 2088, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6987, 0, 3, 1656, 6204, 2133, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 7122, 3, 9, 10, 2181, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7128, 3, 10, 11, 2184, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7134, 3, 11, 12, 2187, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7140, 3, 12, 13, 2190, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7146, 3, 13, 14, 2193, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7152, 3, 14, 15, 2196, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7158, 3, 15, 16, 2199, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7164, 3, 16, 17, 2202, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7170, 3, 17, 18, 2205, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7176, 3, 23, 24, 2211, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7182, 3, 24, 25, 2214, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7188, 3, 25, 26, 2217, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7194, 3, 26, 27, 2220, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7200, 3, 27, 28, 2223, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7206, 3, 28, 29, 2226, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7212, 3, 29, 30, 2229, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7218, 3, 30, 31, 2232, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7224, 3, 31, 32, 2235, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 7230, 0, 3, 2178, 7122, 2247, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7248, 0, 3, 2181, 7128, 2256, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7266, 0, 3, 2184, 7134, 2265, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7284, 0, 3, 2187, 7140, 2274, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7302, 0, 3, 2190, 7146, 2283, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7320, 0, 3, 2193, 7152, 2292, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7338, 0, 3, 2196, 7158, 2301, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7356, 0, 3, 2199, 7164, 2310, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7374, 0, 3, 2202, 7170, 2319, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7392, 0, 3, 2208, 7176, 2337, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7410, 0, 3, 2211, 7182, 2346, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7428, 0, 3, 2214, 7188, 2355, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7446, 0, 3, 2217, 7194, 2364, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7464, 0, 3, 2220, 7200, 2373, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7482, 0, 3, 2223, 7206, 2382, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7500, 0, 3, 2226, 7212, 2391, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7518, 0, 3, 2229, 7218, 2400, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7536, 0, 3, 2232, 7224, 2409, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 7554, 0, 3, 2238, 7230, 100, 106, 2418,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7590, 0, 3, 2247, 7248, 106, 112, 2436,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7626, 0, 3, 2256, 7266, 112, 118, 2454,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7662, 0, 3, 2265, 7284, 118, 124, 2472,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7698, 0, 3, 2274, 7302, 124, 130, 2490,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7734, 0, 3, 2283, 7320, 130, 136, 2508,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7770, 0, 3, 2292, 7338, 136, 142, 2526,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7806, 0, 3, 2301, 7356, 142, 148, 2544,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7842, 0, 3, 2310, 7374, 148, 154, 2562,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7878, 0, 3, 2328, 7392, 166, 172, 2580,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7914, 0, 3, 2337, 7410, 172, 178, 2598,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7950, 0, 3, 2346, 7428, 178, 184, 2616,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7986, 0, 3, 2355, 7446, 184, 190, 2634,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8022, 0, 3, 2364, 7464, 190, 196, 2652,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8058, 0, 3, 2373, 7482, 196, 202, 2670,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8094, 0, 3, 2382, 7500, 202, 208, 2688,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8130, 0, 3, 2391, 7518, 208, 214, 2706,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8166, 0, 3, 2400, 7536, 214, 220, 2724,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8202, 0, 3, 2436, 7626, 232, 242, 2772,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8262, 0, 3, 2454, 7662, 242, 252, 2802,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8322, 0, 3, 2472, 7698, 252, 262, 2832,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8382, 0, 3, 2490, 7734, 262, 272, 2862,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8442, 0, 3, 2508, 7770, 272, 282, 2892,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8502, 0, 3, 2526, 7806, 282, 292, 2922,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8562, 0, 3, 2544, 7842, 292, 302, 2952,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8622, 0, 3, 2598, 7950, 322, 332, 3012,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8682, 0, 3, 2616, 7986, 332, 342, 3042,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8742, 0, 3, 2634, 8022, 342, 352, 3072,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8802, 0, 3, 2652, 8058, 352, 362, 3102,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8862, 0, 3, 2670, 8094, 362, 372, 3132,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8922, 0, 3, 2688, 8130, 372, 382, 3162,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8982, 0, 3, 2706, 8166, 382, 392, 3192,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9042, 0, 3, 7554, 7590, 2742, 8202, 412,
                                                 427, 3222, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9132, 0, 3, 7590, 7626, 2772, 8262, 427,
                                                 442, 3267, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9222, 0, 3, 7626, 7662, 2802, 8322, 442,
                                                 457, 3312, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9312, 0, 3, 7662, 7698, 2832, 8382, 457,
                                                 472, 3357, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9402, 0, 3, 7698, 7734, 2862, 8442, 472,
                                                 487, 3402, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9492, 0, 3, 7734, 7770, 2892, 8502, 487,
                                                 502, 3447, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9582, 0, 3, 7770, 7806, 2922, 8562, 502,
                                                 517, 3492, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9672, 0, 3, 7878, 7914, 2982, 8622, 547,
                                                 562, 3537, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9762, 0, 3, 7914, 7950, 3012, 8682, 562,
                                                 577, 3582, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9852, 0, 3, 7950, 7986, 3042, 8742, 577,
                                                 592, 3627, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9942, 0, 3, 7986, 8022, 3072, 8802, 592,
                                                 607, 3672, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10032, 0, 3, 8022, 8058, 3102, 8862,
                                                 607, 622, 3717, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10122, 0, 3, 8058, 8094, 3132, 8922,
                                                 622, 637, 3762, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10212, 0, 3, 8094, 8130, 3162, 8982,
                                                 637, 652, 3807, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10302, 0, 3, 8202, 8262, 3267, 9222,
                                                 682, 703, 3915, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10428, 0, 3, 8262, 8322, 3312, 9312,
                                                 703, 724, 3978, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10554, 0, 3, 8322, 8382, 3357, 9402,
                                                 724, 745, 4041, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10680, 0, 3, 8382, 8442, 3402, 9492,
                                                 745, 766, 4104, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10806, 0, 3, 8442, 8502, 3447, 9582,
                                                 766, 787, 4167, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10932, 0, 3, 8622, 8682, 3582, 9852,
                                                 829, 850, 4293, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11058, 0, 3, 8682, 8742, 3627, 9942,
                                                 850, 871, 4356, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11184, 0, 3, 8742, 8802, 3672, 10032,
                                                 871, 892, 4419, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11310, 0, 3, 8802, 8862, 3717, 10122,
                                                 892, 913, 4482, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11436, 0, 3, 8862, 8922, 3762, 10212,
                                                 913, 934, 4545, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11562, 0, 3, 9042, 9132, 3852, 10302,
                                                 976, 1004, 4608, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11730, 0, 3, 9132, 9222, 3915, 10428,
                                                 1004, 1032, 4692, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11898, 0, 3, 9222, 9312, 3978, 10554,
                                                 1032, 1060, 4776, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12066, 0, 3, 9312, 9402, 4041, 10680,
                                                 1060, 1088, 4860, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12234, 0, 3, 9402, 9492, 4104, 10806,
                                                 1088, 1116, 4944, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12402, 0, 3, 9672, 9762, 4230, 10932,
                                                 1172, 1200, 5028, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12570, 0, 3, 9762, 9852, 4293, 11058,
                                                 1200, 1228, 5112, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12738, 0, 3, 9852, 9942, 4356, 11184,
                                                 1228, 1256, 5196, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12906, 0, 3, 9942, 10032, 4419, 11310,
                                                 1256, 1284, 5280, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13074, 0, 3, 10032, 10122, 4482, 11436,
                                                 1284, 1312, 5364, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13242, 0, 3, 10302, 10428, 4692, 11898,
                                                 1368, 1404, 5556, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13458, 0, 3, 10428, 10554, 4776, 12066,
                                                 1404, 1440, 5664, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13674, 0, 3, 10554, 10680, 4860, 12234,
                                                 1440, 1476, 5772, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13890, 0, 3, 10932, 11058, 5112, 12738,
                                                 1548, 1584, 5988, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14106, 0, 3, 11058, 11184, 5196, 12906,
                                                 1584, 1620, 6096, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14322, 0, 3, 11184, 11310, 5280, 13074,
                                                 1620, 1656, 6204, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14538, 0, 3, 11562, 11730, 5448, 13242,
                                                 1728, 1773, 6312, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14808, 0, 3, 11730, 11898, 5556, 13458,
                                                 1773, 1818, 6447, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15078, 0, 3, 11898, 12066, 5664, 13674,
                                                 1818, 1863, 6582, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15348, 0, 3, 12402, 12570, 5880, 13890,
                                                 1953, 1998, 6717, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15618, 0, 3, 12570, 12738, 5988, 14106,
                                                 1998, 2043, 6852, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15888, 0, 3, 12738, 12906, 6096, 14322,
                                                 2043, 2088, 6987, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16158, 3, 2178, 2181, 7128, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16168, 3, 2181, 2184, 7134, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16178, 3, 2184, 2187, 7140, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16188, 3, 2187, 2190, 7146, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16198, 3, 2190, 2193, 7152, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16208, 3, 2193, 2196, 7158, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16218, 3, 2196, 2199, 7164, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16228, 3, 2199, 2202, 7170, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16238, 3, 2208, 2211, 7182, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16248, 3, 2211, 2214, 7188, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16258, 3, 2214, 2217, 7194, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16268, 3, 2217, 2220, 7200, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16278, 3, 2220, 2223, 7206, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16288, 3, 2223, 2226, 7212, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16298, 3, 2226, 2229, 7218, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16308, 3, 2229, 2232, 7224, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 16318, 0, 3, 7122, 16158, 7248, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16348, 0, 3, 7128, 16168, 7266, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16378, 0, 3, 7134, 16178, 7284, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16408, 0, 3, 7140, 16188, 7302, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16438, 0, 3, 7146, 16198, 7320, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16468, 0, 3, 7152, 16208, 7338, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16498, 0, 3, 7158, 16218, 7356, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16528, 0, 3, 7164, 16228, 7374, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16558, 0, 3, 7176, 16238, 7410, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16588, 0, 3, 7182, 16248, 7428, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16618, 0, 3, 7188, 16258, 7446, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16648, 0, 3, 7194, 16268, 7464, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16678, 0, 3, 7200, 16278, 7482, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16708, 0, 3, 7206, 16288, 7500, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16738, 0, 3, 7212, 16298, 7518, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16768, 0, 3, 7218, 16308, 7536, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 16798, 0, 3, 7248, 16348, 2418, 2436,
                                                 7626, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16858, 0, 3, 7266, 16378, 2436, 2454,
                                                 7662, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16918, 0, 3, 7284, 16408, 2454, 2472,
                                                 7698, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16978, 0, 3, 7302, 16438, 2472, 2490,
                                                 7734, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17038, 0, 3, 7320, 16468, 2490, 2508,
                                                 7770, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17098, 0, 3, 7338, 16498, 2508, 2526,
                                                 7806, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17158, 0, 3, 7356, 16528, 2526, 2544,
                                                 7842, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17218, 0, 3, 7410, 16588, 2580, 2598,
                                                 7950, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17278, 0, 3, 7428, 16618, 2598, 2616,
                                                 7986, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17338, 0, 3, 7446, 16648, 2616, 2634,
                                                 8022, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17398, 0, 3, 7464, 16678, 2634, 2652,
                                                 8058, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17458, 0, 3, 7482, 16708, 2652, 2670,
                                                 8094, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17518, 0, 3, 7500, 16738, 2670, 2688,
                                                 8130, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17578, 0, 3, 7518, 16768, 2688, 2706,
                                                 8166, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17638, 0, 3, 7626, 16858, 2742, 2772,
                                                 8262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17738, 0, 3, 7662, 16918, 2772, 2802,
                                                 8322, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17838, 0, 3, 7698, 16978, 2802, 2832,
                                                 8382, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17938, 0, 3, 7734, 17038, 2832, 2862,
                                                 8442, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18038, 0, 3, 7770, 17098, 2862, 2892,
                                                 8502, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18138, 0, 3, 7806, 17158, 2892, 2922,
                                                 8562, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18238, 0, 3, 7950, 17278, 2982, 3012,
                                                 8682, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18338, 0, 3, 7986, 17338, 3012, 3042,
                                                 8742, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18438, 0, 3, 8022, 17398, 3042, 3072,
                                                 8802, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18538, 0, 3, 8058, 17458, 3072, 3102,
                                                 8862, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18638, 0, 3, 8094, 17518, 3102, 3132,
                                                 8922, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18738, 0, 3, 8130, 17578, 3132, 3162,
                                                 8982, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18838, 0, 3, 16798, 16858, 8262, 17738,
                                                 3222, 3267, 9222, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18988, 0, 3, 16858, 16918, 8322, 17838,
                                                 3267, 3312, 9312, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19138, 0, 3, 16918, 16978, 8382, 17938,
                                                 3312, 3357, 9402, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19288, 0, 3, 16978, 17038, 8442, 18038,
                                                 3357, 3402, 9492, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19438, 0, 3, 17038, 17098, 8502, 18138,
                                                 3402, 3447, 9582, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19588, 0, 3, 17218, 17278, 8682, 18338,
                                                 3537, 3582, 9852, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19738, 0, 3, 17278, 17338, 8742, 18438,
                                                 3582, 3627, 9942, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19888, 0, 3, 17338, 17398, 8802, 18538,
                                                 3627, 3672, 10032, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20038, 0, 3, 17398, 17458, 8862, 18638,
                                                 3672, 3717, 10122, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20188, 0, 3, 17458, 17518, 8922, 18738,
                                                 3717, 3762, 10212, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20338, 0, 3, 17638, 17738, 9222, 18988,
                                                 3852, 3915, 10428, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20548, 0, 3, 17738, 17838, 9312, 19138,
                                                 3915, 3978, 10554, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20758, 0, 3, 17838, 17938, 9402, 19288,
                                                 3978, 4041, 10680, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20968, 0, 3, 17938, 18038, 9492, 19438,
                                                 4041, 4104, 10806, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21178, 0, 3, 18238, 18338, 9852, 19738,
                                                 4230, 4293, 11058, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21388, 0, 3, 18338, 18438, 9942, 19888,
                                                 4293, 4356, 11184, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21598, 0, 3, 18438, 18538, 10032, 20038,
                                                 4356, 4419, 11310, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21808, 0, 3, 18538, 18638, 10122, 20188,
                                                 4419, 4482, 11436, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22018, 0, 3, 18838, 18988, 10428, 20548,
                                                 4608, 4692, 11898, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22298, 0, 3, 18988, 19138, 10554, 20758,
                                                 4692, 4776, 12066, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22578, 0, 3, 19138, 19288, 10680, 20968,
                                                 4776, 4860, 12234, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22858, 0, 3, 19588, 19738, 11058, 21388,
                                                 5028, 5112, 12738, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23138, 0, 3, 19738, 19888, 11184, 21598,
                                                 5112, 5196, 12906, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23418, 0, 3, 19888, 20038, 11310, 21808,
                                                 5196, 5280, 13074, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 23698, 0, 3, 20338, 20548, 11898, 22298,
                                                 5448, 5556, 13458, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 24058, 0, 3, 20548, 20758, 12066, 22578,
                                                 5556, 5664, 13674, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 24418, 0, 3, 21178, 21388, 12738, 23138,
                                                 5880, 5988, 14106, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 24778, 0, 3, 21388, 21598, 12906, 23418,
                                                 5988, 6096, 14322, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 25138, 0, 3, 22018, 22298, 13458, 24058,
                                                 6312, 6447, 15078, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 25588, 0, 3, 22858, 23138, 14106, 24778,
                                                 6717, 6852, 15888, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26038, 3, 7122, 7128, 16168, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26053, 3, 7128, 7134, 16178, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26068, 3, 7134, 7140, 16188, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26083, 3, 7140, 7146, 16198, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26098, 3, 7146, 7152, 16208, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26113, 3, 7152, 7158, 16218, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26128, 3, 7158, 7164, 16228, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26143, 3, 7176, 7182, 16248, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26158, 3, 7182, 7188, 16258, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26173, 3, 7188, 7194, 16268, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26188, 3, 7194, 7200, 16278, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26203, 3, 7200, 7206, 16288, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26218, 3, 7206, 7212, 16298, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26233, 3, 7212, 7218, 16308, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 26248, 0, 3, 16158, 26038, 16348, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26293, 0, 3, 16168, 26053, 16378, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26338, 0, 3, 16178, 26068, 16408, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26383, 0, 3, 16188, 26083, 16438, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26428, 0, 3, 16198, 26098, 16468, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26473, 0, 3, 16208, 26113, 16498, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26518, 0, 3, 16218, 26128, 16528, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26563, 0, 3, 16238, 26143, 16588, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26608, 0, 3, 16248, 26158, 16618, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26653, 0, 3, 16258, 26173, 16648, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26698, 0, 3, 16268, 26188, 16678, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26743, 0, 3, 16278, 26203, 16708, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26788, 0, 3, 16288, 26218, 16738, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26833, 0, 3, 16298, 26233, 16768, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 26878, 0, 3, 16318, 26248, 7554, 7590,
                                                 16798, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26968, 0, 3, 16348, 26293, 7590, 7626,
                                                 16858, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27058, 0, 3, 16378, 26338, 7626, 7662,
                                                 16918, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27148, 0, 3, 16408, 26383, 7662, 7698,
                                                 16978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27238, 0, 3, 16438, 26428, 7698, 7734,
                                                 17038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27328, 0, 3, 16468, 26473, 7734, 7770,
                                                 17098, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27418, 0, 3, 16498, 26518, 7770, 7806,
                                                 17158, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27508, 0, 3, 16558, 26563, 7878, 7914,
                                                 17218, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27598, 0, 3, 16588, 26608, 7914, 7950,
                                                 17278, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27688, 0, 3, 16618, 26653, 7950, 7986,
                                                 17338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27778, 0, 3, 16648, 26698, 7986, 8022,
                                                 17398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27868, 0, 3, 16678, 26743, 8022, 8058,
                                                 17458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27958, 0, 3, 16708, 26788, 8058, 8094,
                                                 17518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28048, 0, 3, 16738, 26833, 8094, 8130,
                                                 17578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28138, 0, 3, 16858, 27058, 8202, 8262,
                                                 17738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28288, 0, 3, 16918, 27148, 8262, 8322,
                                                 17838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28438, 0, 3, 16978, 27238, 8322, 8382,
                                                 17938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28588, 0, 3, 17038, 27328, 8382, 8442,
                                                 18038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28738, 0, 3, 17098, 27418, 8442, 8502,
                                                 18138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28888, 0, 3, 17278, 27688, 8622, 8682,
                                                 18338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29038, 0, 3, 17338, 27778, 8682, 8742,
                                                 18438, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29188, 0, 3, 17398, 27868, 8742, 8802,
                                                 18538, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29338, 0, 3, 17458, 27958, 8802, 8862,
                                                 18638, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29488, 0, 3, 17518, 28048, 8862, 8922,
                                                 18738, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29638, 0, 3, 26878, 26968, 17638, 28138,
                                                 9042, 9132, 18838, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29863, 0, 3, 26968, 27058, 17738, 28288,
                                                 9132, 9222, 18988, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30088, 0, 3, 27058, 27148, 17838, 28438,
                                                 9222, 9312, 19138, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30313, 0, 3, 27148, 27238, 17938, 28588,
                                                 9312, 9402, 19288, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30538, 0, 3, 27238, 27328, 18038, 28738,
                                                 9402, 9492, 19438, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30763, 0, 3, 27508, 27598, 18238, 28888,
                                                 9672, 9762, 19588, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30988, 0, 3, 27598, 27688, 18338, 29038,
                                                 9762, 9852, 19738, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31213, 0, 3, 27688, 27778, 18438, 29188,
                                                 9852, 9942, 19888, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31438, 0, 3, 27778, 27868, 18538, 29338,
                                                 9942, 10032, 20038, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31663, 0, 3, 27868, 27958, 18638, 29488,
                                                 10032, 10122, 20188, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31888, 0, 3, 28138, 28288, 18988, 30088,
                                                 10302, 10428, 20548, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32203, 0, 3, 28288, 28438, 19138, 30313,
                                                 10428, 10554, 20758, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32518, 0, 3, 28438, 28588, 19288, 30538,
                                                 10554, 10680, 20968, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32833, 0, 3, 28888, 29038, 19738, 31213,
                                                 10932, 11058, 21388, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33148, 0, 3, 29038, 29188, 19888, 31438,
                                                 11058, 11184, 21598, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33463, 0, 3, 29188, 29338, 20038, 31663,
                                                 11184, 11310, 21808, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 33778, 0, 3, 29638, 29863, 20338, 31888,
                                                 11562, 11730, 22018, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34198, 0, 3, 29863, 30088, 20548, 32203,
                                                 11730, 11898, 22298, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34618, 0, 3, 30088, 30313, 20758, 32518,
                                                 11898, 12066, 22578, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35038, 0, 3, 30763, 30988, 21178, 32833,
                                                 12402, 12570, 22858, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35458, 0, 3, 30988, 31213, 21388, 33148,
                                                 12570, 12738, 23138, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35878, 0, 3, 31213, 31438, 21598, 33463,
                                                 12738, 12906, 23418, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36298, 0, 3, 31888, 32203, 22298, 34618,
                                                 13242, 13458, 24058, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36838, 0, 3, 32833, 33148, 23138, 35878,
                                                 13890, 14106, 24778, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 37378, 0, 3, 33778, 34198, 23698, 36298,
                                                 14538, 14808, 25138, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 38053, 0, 3, 35038, 35458, 24418, 36838,
                                                 15348, 15618, 25588, ncols, alpha, beta, p);

            simdgeo::geom_k_x(buffer, 38728, 35038, 38053, 1, 15, ncols, alpha);

            simdgeo::geom_k_y(buffer, 39268, 35038, 38053, 1, 15, ncols, alpha);

            simdgeo::geom_k_z(buffer, 39808, 35038, 38053, 1, 15, ncols, alpha);

            simdgeo::geom_k_x(buffer, 40348, 33778, 37378, 1, 15, ncols, alpha);

            simdgeo::geom_k_y(buffer, 40888, 33778, 37378, 1, 15, ncols, alpha);

            simdgeo::geom_k_z(buffer, 41428, 33778, 37378, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 41968, 40348, 1620, ncols);

            simdfunc::contract_primitives(buffer, 43588, 38728, 1620, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 45208, 43588, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 45208, 9, nmax);

    simdtrf::transform_g_inner(buffer, 45208, 44128, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 135 * nvalues, nvalues, buffer, 45208, 9, nmax);

    simdtrf::transform_g_inner(buffer, 45208, 44668, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 270 * nvalues, nvalues, buffer, 45208, 9, nmax);

    simdtrf::transform_g_inner(buffer, 45208, 41968, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 405 * nvalues, nvalues, buffer, 45208, 9, nmax);

    simdtrf::transform_g_inner(buffer, 45208, 42508, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 540 * nvalues, nvalues, buffer, 45208, 9, nmax);

    simdtrf::transform_g_inner(buffer, 45208, 43048, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 675 * nvalues, nvalues, buffer, 45208, 9, nmax);
}

}  // namespace simdt2ceri
