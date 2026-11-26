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


#include "SimdElectronRepulsionGeom10RsRecLF.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_lf_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_lf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 35105, 32090, 2700, nvalues);

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

            compute_prim_is_electron_repulsion_0(buffer, 932, 0, 398, 413, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 960, 0, 413, 428, 680, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 988, 0, 428, 443, 701, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1016, 0, 443, 458, 722, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1044, 0, 458, 473, 743, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1072, 0, 473, 488, 764, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1100, 0, 518, 533, 806, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1128, 0, 533, 548, 827, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1156, 0, 548, 563, 848, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1184, 0, 563, 578, 869, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1212, 0, 578, 593, 890, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1240, 0, 593, 608, 911, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1268, 0, 638, 659, 960, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1304, 0, 659, 680, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1340, 0, 680, 701, 1016, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1376, 0, 701, 722, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1412, 0, 722, 743, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1448, 0, 785, 806, 1128, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1484, 0, 806, 827, 1156, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1520, 0, 827, 848, 1184, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1556, 0, 848, 869, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1592, 0, 869, 890, 1240, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1628, 0, 932, 960, 1304, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1673, 0, 960, 988, 1340, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1718, 0, 988, 1016, 1376, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1763, 0, 1016, 1044, 1412, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1808, 0, 1100, 1128, 1484, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1853, 0, 1128, 1156, 1520, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1898, 0, 1156, 1184, 1556, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1943, 0, 1184, 1212, 1592, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1988, 0, 1268, 1304, 1673, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2043, 0, 1304, 1340, 1718, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2098, 0, 1340, 1376, 1763, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2153, 0, 1448, 1484, 1853, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2208, 0, 1484, 1520, 1898, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2263, 0, 1520, 1556, 1943, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2318, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2321, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2324, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2327, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2330, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2333, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2336, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2339, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2342, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2345, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2348, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2351, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2354, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2357, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2360, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2363, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2366, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2369, 3, 31, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2372, 3, 9, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2381, 3, 10, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2390, 3, 11, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2399, 3, 12, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2408, 3, 13, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2417, 3, 14, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2426, 3, 15, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2435, 3, 16, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2444, 3, 17, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2453, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2462, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2471, 3, 24, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2480, 3, 25, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2489, 3, 26, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2498, 3, 27, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2507, 3, 28, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2516, 3, 29, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2525, 3, 30, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2534, 0, 3, 35, 2372, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2552, 0, 3, 38, 2381, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2570, 0, 3, 41, 2390, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2588, 0, 3, 44, 2399, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2606, 0, 3, 47, 2408, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2624, 0, 3, 50, 2417, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2642, 0, 3, 53, 2426, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2660, 0, 3, 56, 2435, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2678, 0, 3, 59, 2444, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2696, 0, 3, 68, 2453, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2714, 0, 3, 71, 2462, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2732, 0, 3, 74, 2471, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2750, 0, 3, 77, 2480, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2768, 0, 3, 80, 2489, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2786, 0, 3, 83, 2498, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2804, 0, 3, 86, 2507, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2822, 0, 3, 89, 2516, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2840, 0, 3, 92, 2525, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2858, 0, 3, 98, 2534, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2888, 0, 3, 104, 2552, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2918, 0, 3, 110, 2570, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2948, 0, 3, 116, 2588, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2978, 0, 3, 122, 2606, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3008, 0, 3, 128, 2624, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3038, 0, 3, 134, 2642, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3068, 0, 3, 140, 2660, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3098, 0, 3, 146, 2678, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3128, 0, 3, 158, 2696, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3158, 0, 3, 164, 2714, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3188, 0, 3, 170, 2732, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3218, 0, 3, 176, 2750, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3248, 0, 3, 182, 2768, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3278, 0, 3, 188, 2786, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3308, 0, 3, 194, 2804, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3338, 0, 3, 200, 2822, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3368, 0, 3, 206, 2840, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3398, 0, 3, 228, 2918, 413, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3443, 0, 3, 238, 2948, 428, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3488, 0, 3, 248, 2978, 443, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3533, 0, 3, 258, 3008, 458, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3578, 0, 3, 268, 3038, 473, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3623, 0, 3, 278, 3068, 488, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3668, 0, 3, 288, 3098, 503, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3713, 0, 3, 318, 3188, 533, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3758, 0, 3, 328, 3218, 548, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3803, 0, 3, 338, 3248, 563, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3848, 0, 3, 348, 3278, 578, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3893, 0, 3, 358, 3308, 593, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3938, 0, 3, 368, 3338, 608, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3983, 0, 3, 378, 3368, 623, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4028, 0, 3, 398, 3398, 638, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4091, 0, 3, 413, 3443, 659, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4154, 0, 3, 428, 3488, 680, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4217, 0, 3, 443, 3533, 701, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4280, 0, 3, 458, 3578, 722, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4343, 0, 3, 473, 3623, 743, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4406, 0, 3, 488, 3668, 764, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4469, 0, 3, 518, 3713, 785, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4532, 0, 3, 533, 3758, 806, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4595, 0, 3, 548, 3803, 827, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4658, 0, 3, 563, 3848, 848, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4721, 0, 3, 578, 3893, 869, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4784, 0, 3, 593, 3938, 890, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4847, 0, 3, 608, 3983, 911, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4910, 0, 3, 659, 4154, 960, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4994, 0, 3, 680, 4217, 988, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5078, 0, 3, 701, 4280, 1016, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5162, 0, 3, 722, 4343, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5246, 0, 3, 743, 4406, 1072, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5330, 0, 3, 806, 4595, 1128, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5414, 0, 3, 827, 4658, 1156, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5498, 0, 3, 848, 4721, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5582, 0, 3, 869, 4784, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5666, 0, 3, 890, 4847, 1240, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5750, 0, 3, 932, 4910, 1268, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5858, 0, 3, 960, 4994, 1304, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5966, 0, 3, 988, 5078, 1340, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 6074, 0, 3, 1016, 5162, 1376, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6182, 0, 3, 1044, 5246, 1412, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6290, 0, 3, 1100, 5330, 1448, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6398, 0, 3, 1128, 5414, 1484, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6506, 0, 3, 1156, 5498, 1520, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6614, 0, 3, 1184, 5582, 1556, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6722, 0, 3, 1212, 5666, 1592, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6830, 0, 3, 1304, 5966, 1673, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6965, 0, 3, 1340, 6074, 1718, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7100, 0, 3, 1376, 6182, 1763, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7235, 0, 3, 1484, 6506, 1853, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7370, 0, 3, 1520, 6614, 1898, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7505, 0, 3, 1556, 6722, 1943, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7640, 0, 3, 1628, 6830, 1988, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7805, 0, 3, 1673, 6965, 2043, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7970, 0, 3, 1718, 7100, 2098, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8135, 0, 3, 1808, 7235, 2153, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8300, 0, 3, 1853, 7370, 2208, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8465, 0, 3, 1898, 7505, 2263, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8630, 3, 9, 10, 2321, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8636, 3, 10, 11, 2324, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8642, 3, 11, 12, 2327, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8648, 3, 12, 13, 2330, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8654, 3, 13, 14, 2333, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8660, 3, 14, 15, 2336, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8666, 3, 15, 16, 2339, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8672, 3, 16, 17, 2342, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8678, 3, 22, 23, 2348, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8684, 3, 23, 24, 2351, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8690, 3, 24, 25, 2354, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8696, 3, 25, 26, 2357, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8702, 3, 26, 27, 2360, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8708, 3, 27, 28, 2363, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8714, 3, 28, 29, 2366, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8720, 3, 29, 30, 2369, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 8726, 0, 3, 2318, 8630, 2381, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8744, 0, 3, 2321, 8636, 2390, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8762, 0, 3, 2324, 8642, 2399, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8780, 0, 3, 2327, 8648, 2408, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8798, 0, 3, 2330, 8654, 2417, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8816, 0, 3, 2333, 8660, 2426, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8834, 0, 3, 2336, 8666, 2435, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8852, 0, 3, 2339, 8672, 2444, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8870, 0, 3, 2345, 8678, 2462, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8888, 0, 3, 2348, 8684, 2471, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8906, 0, 3, 2351, 8690, 2480, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8924, 0, 3, 2354, 8696, 2489, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8942, 0, 3, 2357, 8702, 2498, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8960, 0, 3, 2360, 8708, 2507, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8978, 0, 3, 2363, 8714, 2516, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8996, 0, 3, 2366, 8720, 2525, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 9014, 0, 3, 2372, 8726, 98, 104, 2552,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9050, 0, 3, 2381, 8744, 104, 110, 2570,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9086, 0, 3, 2390, 8762, 110, 116, 2588,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9122, 0, 3, 2399, 8780, 116, 122, 2606,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9158, 0, 3, 2408, 8798, 122, 128, 2624,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9194, 0, 3, 2417, 8816, 128, 134, 2642,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9230, 0, 3, 2426, 8834, 134, 140, 2660,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9266, 0, 3, 2435, 8852, 140, 146, 2678,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9302, 0, 3, 2453, 8870, 158, 164, 2714,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9338, 0, 3, 2462, 8888, 164, 170, 2732,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9374, 0, 3, 2471, 8906, 170, 176, 2750,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9410, 0, 3, 2480, 8924, 176, 182, 2768,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9446, 0, 3, 2489, 8942, 182, 188, 2786,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9482, 0, 3, 2498, 8960, 188, 194, 2804,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9518, 0, 3, 2507, 8978, 194, 200, 2822,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9554, 0, 3, 2516, 8996, 200, 206, 2840,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9590, 0, 3, 2552, 9050, 218, 228, 2918,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9650, 0, 3, 2570, 9086, 228, 238, 2948,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9710, 0, 3, 2588, 9122, 238, 248, 2978,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9770, 0, 3, 2606, 9158, 248, 258, 3008,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9830, 0, 3, 2624, 9194, 258, 268, 3038,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9890, 0, 3, 2642, 9230, 268, 278, 3068,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9950, 0, 3, 2660, 9266, 278, 288, 3098,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10010, 0, 3, 2714, 9338, 308, 318, 3188,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10070, 0, 3, 2732, 9374, 318, 328, 3218,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10130, 0, 3, 2750, 9410, 328, 338, 3248,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10190, 0, 3, 2768, 9446, 338, 348, 3278,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10250, 0, 3, 2786, 9482, 348, 358, 3308,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10310, 0, 3, 2804, 9518, 358, 368, 3338,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10370, 0, 3, 2822, 9554, 368, 378, 3368,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10430, 0, 3, 9014, 9050, 2918, 9650,
                                                 398, 413, 3443, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10520, 0, 3, 9050, 9086, 2948, 9710,
                                                 413, 428, 3488, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10610, 0, 3, 9086, 9122, 2978, 9770,
                                                 428, 443, 3533, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10700, 0, 3, 9122, 9158, 3008, 9830,
                                                 443, 458, 3578, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10790, 0, 3, 9158, 9194, 3038, 9890,
                                                 458, 473, 3623, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10880, 0, 3, 9194, 9230, 3068, 9950,
                                                 473, 488, 3668, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10970, 0, 3, 9302, 9338, 3188, 10070,
                                                 518, 533, 3758, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11060, 0, 3, 9338, 9374, 3218, 10130,
                                                 533, 548, 3803, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11150, 0, 3, 9374, 9410, 3248, 10190,
                                                 548, 563, 3848, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11240, 0, 3, 9410, 9446, 3278, 10250,
                                                 563, 578, 3893, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11330, 0, 3, 9446, 9482, 3308, 10310,
                                                 578, 593, 3938, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11420, 0, 3, 9482, 9518, 3338, 10370,
                                                 593, 608, 3983, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11510, 0, 3, 9590, 9650, 3443, 10520,
                                                 638, 659, 4154, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11636, 0, 3, 9650, 9710, 3488, 10610,
                                                 659, 680, 4217, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11762, 0, 3, 9710, 9770, 3533, 10700,
                                                 680, 701, 4280, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11888, 0, 3, 9770, 9830, 3578, 10790,
                                                 701, 722, 4343, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12014, 0, 3, 9830, 9890, 3623, 10880,
                                                 722, 743, 4406, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12140, 0, 3, 10010, 10070, 3758, 11060,
                                                 785, 806, 4595, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12266, 0, 3, 10070, 10130, 3803, 11150,
                                                 806, 827, 4658, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12392, 0, 3, 10130, 10190, 3848, 11240,
                                                 827, 848, 4721, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12518, 0, 3, 10190, 10250, 3893, 11330,
                                                 848, 869, 4784, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12644, 0, 3, 10250, 10310, 3938, 11420,
                                                 869, 890, 4847, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12770, 0, 3, 10430, 10520, 4154, 11636,
                                                 932, 960, 4994, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12938, 0, 3, 10520, 10610, 4217, 11762,
                                                 960, 988, 5078, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13106, 0, 3, 10610, 10700, 4280, 11888,
                                                 988, 1016, 5162, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13274, 0, 3, 10700, 10790, 4343, 12014,
                                                 1016, 1044, 5246, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13442, 0, 3, 10970, 11060, 4595, 12266,
                                                 1100, 1128, 5414, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13610, 0, 3, 11060, 11150, 4658, 12392,
                                                 1128, 1156, 5498, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13778, 0, 3, 11150, 11240, 4721, 12518,
                                                 1156, 1184, 5582, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13946, 0, 3, 11240, 11330, 4784, 12644,
                                                 1184, 1212, 5666, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14114, 0, 3, 11510, 11636, 4994, 12938,
                                                 1268, 1304, 5966, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14330, 0, 3, 11636, 11762, 5078, 13106,
                                                 1304, 1340, 6074, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14546, 0, 3, 11762, 11888, 5162, 13274,
                                                 1340, 1376, 6182, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14762, 0, 3, 12140, 12266, 5414, 13610,
                                                 1448, 1484, 6506, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14978, 0, 3, 12266, 12392, 5498, 13778,
                                                 1484, 1520, 6614, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15194, 0, 3, 12392, 12518, 5582, 13946,
                                                 1520, 1556, 6722, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15410, 0, 3, 12770, 12938, 5966, 14330,
                                                 1628, 1673, 6965, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15680, 0, 3, 12938, 13106, 6074, 14546,
                                                 1673, 1718, 7100, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15950, 0, 3, 13442, 13610, 6506, 14978,
                                                 1808, 1853, 7370, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 16220, 0, 3, 13610, 13778, 6614, 15194,
                                                 1853, 1898, 7505, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 16490, 0, 3, 14114, 14330, 6965, 15680,
                                                 1988, 2043, 7970, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 16820, 0, 3, 14762, 14978, 7370, 16220,
                                                 2153, 2208, 8465, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17150, 3, 2318, 2321, 8636, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17160, 3, 2321, 2324, 8642, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17170, 3, 2324, 2327, 8648, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17180, 3, 2327, 2330, 8654, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17190, 3, 2330, 2333, 8660, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17200, 3, 2333, 2336, 8666, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17210, 3, 2336, 2339, 8672, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17220, 3, 2345, 2348, 8684, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17230, 3, 2348, 2351, 8690, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17240, 3, 2351, 2354, 8696, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17250, 3, 2354, 2357, 8702, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17260, 3, 2357, 2360, 8708, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17270, 3, 2360, 2363, 8714, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17280, 3, 2363, 2366, 8720, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 17290, 0, 3, 8630, 17150, 8744, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17320, 0, 3, 8636, 17160, 8762, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17350, 0, 3, 8642, 17170, 8780, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17380, 0, 3, 8648, 17180, 8798, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17410, 0, 3, 8654, 17190, 8816, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17440, 0, 3, 8660, 17200, 8834, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17470, 0, 3, 8666, 17210, 8852, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17500, 0, 3, 8678, 17220, 8888, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17530, 0, 3, 8684, 17230, 8906, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17560, 0, 3, 8690, 17240, 8924, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17590, 0, 3, 8696, 17250, 8942, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17620, 0, 3, 8702, 17260, 8960, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17650, 0, 3, 8708, 17270, 8978, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17680, 0, 3, 8714, 17280, 8996, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 17710, 0, 3, 8726, 17290, 2534, 2552,
                                                 9050, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17770, 0, 3, 8744, 17320, 2552, 2570,
                                                 9086, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17830, 0, 3, 8762, 17350, 2570, 2588,
                                                 9122, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17890, 0, 3, 8780, 17380, 2588, 2606,
                                                 9158, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17950, 0, 3, 8798, 17410, 2606, 2624,
                                                 9194, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18010, 0, 3, 8816, 17440, 2624, 2642,
                                                 9230, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18070, 0, 3, 8834, 17470, 2642, 2660,
                                                 9266, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18130, 0, 3, 8870, 17500, 2696, 2714,
                                                 9338, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18190, 0, 3, 8888, 17530, 2714, 2732,
                                                 9374, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18250, 0, 3, 8906, 17560, 2732, 2750,
                                                 9410, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18310, 0, 3, 8924, 17590, 2750, 2768,
                                                 9446, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18370, 0, 3, 8942, 17620, 2768, 2786,
                                                 9482, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18430, 0, 3, 8960, 17650, 2786, 2804,
                                                 9518, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18490, 0, 3, 8978, 17680, 2804, 2822,
                                                 9554, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18550, 0, 3, 9014, 17710, 2858, 2888,
                                                 9590, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18650, 0, 3, 9050, 17770, 2888, 2918,
                                                 9650, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18750, 0, 3, 9086, 17830, 2918, 2948,
                                                 9710, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18850, 0, 3, 9122, 17890, 2948, 2978,
                                                 9770, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18950, 0, 3, 9158, 17950, 2978, 3008,
                                                 9830, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19050, 0, 3, 9194, 18010, 3008, 3038,
                                                 9890, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19150, 0, 3, 9230, 18070, 3038, 3068,
                                                 9950, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19250, 0, 3, 9302, 18130, 3128, 3158,
                                                 10010, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19350, 0, 3, 9338, 18190, 3158, 3188,
                                                 10070, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19450, 0, 3, 9374, 18250, 3188, 3218,
                                                 10130, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19550, 0, 3, 9410, 18310, 3218, 3248,
                                                 10190, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19650, 0, 3, 9446, 18370, 3248, 3278,
                                                 10250, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19750, 0, 3, 9482, 18430, 3278, 3308,
                                                 10310, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19850, 0, 3, 9518, 18490, 3308, 3338,
                                                 10370, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19950, 0, 3, 17710, 17770, 9650, 18750,
                                                 3398, 3443, 10520, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20100, 0, 3, 17770, 17830, 9710, 18850,
                                                 3443, 3488, 10610, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20250, 0, 3, 17830, 17890, 9770, 18950,
                                                 3488, 3533, 10700, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20400, 0, 3, 17890, 17950, 9830, 19050,
                                                 3533, 3578, 10790, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20550, 0, 3, 17950, 18010, 9890, 19150,
                                                 3578, 3623, 10880, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20700, 0, 3, 18130, 18190, 10070, 19450,
                                                 3713, 3758, 11060, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20850, 0, 3, 18190, 18250, 10130, 19550,
                                                 3758, 3803, 11150, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21000, 0, 3, 18250, 18310, 10190, 19650,
                                                 3803, 3848, 11240, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21150, 0, 3, 18310, 18370, 10250, 19750,
                                                 3848, 3893, 11330, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21300, 0, 3, 18370, 18430, 10310, 19850,
                                                 3893, 3938, 11420, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21450, 0, 3, 18550, 18650, 10430, 19950,
                                                 4028, 4091, 11510, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21660, 0, 3, 18650, 18750, 10520, 20100,
                                                 4091, 4154, 11636, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21870, 0, 3, 18750, 18850, 10610, 20250,
                                                 4154, 4217, 11762, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22080, 0, 3, 18850, 18950, 10700, 20400,
                                                 4217, 4280, 11888, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22290, 0, 3, 18950, 19050, 10790, 20550,
                                                 4280, 4343, 12014, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22500, 0, 3, 19250, 19350, 10970, 20700,
                                                 4469, 4532, 12140, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22710, 0, 3, 19350, 19450, 11060, 20850,
                                                 4532, 4595, 12266, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22920, 0, 3, 19450, 19550, 11150, 21000,
                                                 4595, 4658, 12392, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23130, 0, 3, 19550, 19650, 11240, 21150,
                                                 4658, 4721, 12518, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23340, 0, 3, 19650, 19750, 11330, 21300,
                                                 4721, 4784, 12644, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23550, 0, 3, 19950, 20100, 11636, 21870,
                                                 4910, 4994, 12938, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23830, 0, 3, 20100, 20250, 11762, 22080,
                                                 4994, 5078, 13106, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24110, 0, 3, 20250, 20400, 11888, 22290,
                                                 5078, 5162, 13274, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24390, 0, 3, 20700, 20850, 12266, 22920,
                                                 5330, 5414, 13610, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24670, 0, 3, 20850, 21000, 12392, 23130,
                                                 5414, 5498, 13778, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24950, 0, 3, 21000, 21150, 12518, 23340,
                                                 5498, 5582, 13946, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25230, 0, 3, 21450, 21660, 12770, 23550,
                                                 5750, 5858, 14114, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25590, 0, 3, 21660, 21870, 12938, 23830,
                                                 5858, 5966, 14330, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25950, 0, 3, 21870, 22080, 13106, 24110,
                                                 5966, 6074, 14546, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26310, 0, 3, 22500, 22710, 13442, 24390,
                                                 6290, 6398, 14762, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26670, 0, 3, 22710, 22920, 13610, 24670,
                                                 6398, 6506, 14978, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 27030, 0, 3, 22920, 23130, 13778, 24950,
                                                 6506, 6614, 15194, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 27390, 0, 3, 23550, 23830, 14330, 25950,
                                                 6830, 6965, 15680, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 27840, 0, 3, 24390, 24670, 14978, 27030,
                                                 7235, 7370, 16220, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 28290, 0, 3, 25230, 25590, 15410, 27390,
                                                 7640, 7805, 16490, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 28840, 0, 3, 26310, 26670, 15950, 27840,
                                                 8135, 8300, 16820, ncols, alpha, beta, p);

            simdgeo::geom_l_x(buffer, 29390, 26310, 28840, 1, 10, ncols, alpha);

            simdgeo::geom_l_y(buffer, 29840, 26310, 28840, 1, 10, ncols, alpha);

            simdgeo::geom_l_z(buffer, 30290, 26310, 28840, 1, 10, ncols, alpha);

            simdgeo::geom_l_x(buffer, 30740, 25230, 28290, 1, 10, ncols, alpha);

            simdgeo::geom_l_y(buffer, 31190, 25230, 28290, 1, 10, ncols, alpha);

            simdgeo::geom_l_z(buffer, 31640, 25230, 28290, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 32090, 30740, 1350, ncols);

            simdfunc::contract_primitives(buffer, 33440, 29390, 1350, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 34790, 33440, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 34790, 7, nmax);

    simdtrf::transform_f_inner(buffer, 34790, 33890, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 119 * nvalues, nvalues, buffer, 34790, 7, nmax);

    simdtrf::transform_f_inner(buffer, 34790, 34340, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 238 * nvalues, nvalues, buffer, 34790, 7, nmax);

    simdtrf::transform_f_inner(buffer, 34790, 32090, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 357 * nvalues, nvalues, buffer, 34790, 7, nmax);

    simdtrf::transform_f_inner(buffer, 34790, 32540, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 476 * nvalues, nvalues, buffer, 34790, 7, nmax);

    simdtrf::transform_f_inner(buffer, 34790, 32990, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 595 * nvalues, nvalues, buffer, 34790, 7, nmax);
}

}  // namespace simdt2ceri
