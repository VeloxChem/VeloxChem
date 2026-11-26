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


#include "SimdElectronRepulsionGeom10RsRecIH.hpp"

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
#include "SimdGeometryI1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ih_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ih_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 51698, 47862, 3528, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 1628, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1631, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1634, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1637, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1640, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1643, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1646, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1649, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1652, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1655, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1658, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1661, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1664, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1667, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1670, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1673, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1676, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1679, 3, 31, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1682, 3, 9, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1691, 3, 10, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1700, 3, 11, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1709, 3, 12, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1718, 3, 13, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1727, 3, 14, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1736, 3, 15, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1745, 3, 16, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1754, 3, 17, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1763, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1772, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1781, 3, 24, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1790, 3, 25, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1799, 3, 26, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1808, 3, 27, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1817, 3, 28, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1826, 3, 29, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1835, 3, 30, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1844, 0, 3, 35, 1682, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1862, 0, 3, 38, 1691, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1880, 0, 3, 41, 1700, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1898, 0, 3, 44, 1709, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1916, 0, 3, 47, 1718, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1934, 0, 3, 50, 1727, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1952, 0, 3, 53, 1736, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1970, 0, 3, 56, 1745, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1988, 0, 3, 59, 1754, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2006, 0, 3, 68, 1763, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2024, 0, 3, 71, 1772, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2042, 0, 3, 74, 1781, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2060, 0, 3, 77, 1790, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2078, 0, 3, 80, 1799, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2096, 0, 3, 83, 1808, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2114, 0, 3, 86, 1817, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2132, 0, 3, 89, 1826, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2150, 0, 3, 92, 1835, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2168, 0, 3, 98, 1844, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2198, 0, 3, 104, 1862, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2228, 0, 3, 110, 1880, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2258, 0, 3, 116, 1898, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2288, 0, 3, 122, 1916, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2318, 0, 3, 128, 1934, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2348, 0, 3, 134, 1952, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2378, 0, 3, 140, 1970, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2408, 0, 3, 146, 1988, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2438, 0, 3, 158, 2006, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2468, 0, 3, 164, 2024, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2498, 0, 3, 170, 2042, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2528, 0, 3, 176, 2060, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2558, 0, 3, 182, 2078, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2588, 0, 3, 188, 2096, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2618, 0, 3, 194, 2114, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2648, 0, 3, 200, 2132, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2678, 0, 3, 206, 2150, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2708, 0, 3, 228, 2228, 413, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2753, 0, 3, 238, 2258, 428, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2798, 0, 3, 248, 2288, 443, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2843, 0, 3, 258, 2318, 458, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2888, 0, 3, 268, 2348, 473, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2933, 0, 3, 278, 2378, 488, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2978, 0, 3, 288, 2408, 503, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3023, 0, 3, 318, 2498, 533, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3068, 0, 3, 328, 2528, 548, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3113, 0, 3, 338, 2558, 563, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3158, 0, 3, 348, 2588, 578, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3203, 0, 3, 358, 2618, 593, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3248, 0, 3, 368, 2648, 608, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3293, 0, 3, 378, 2678, 623, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3338, 0, 3, 398, 2708, 638, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3401, 0, 3, 413, 2753, 659, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3464, 0, 3, 428, 2798, 680, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3527, 0, 3, 443, 2843, 701, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3590, 0, 3, 458, 2888, 722, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3653, 0, 3, 473, 2933, 743, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3716, 0, 3, 488, 2978, 764, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3779, 0, 3, 518, 3023, 785, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3842, 0, 3, 533, 3068, 806, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3905, 0, 3, 548, 3113, 827, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3968, 0, 3, 563, 3158, 848, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4031, 0, 3, 578, 3203, 869, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4094, 0, 3, 593, 3248, 890, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4157, 0, 3, 608, 3293, 911, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4220, 0, 3, 659, 3464, 960, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4304, 0, 3, 680, 3527, 988, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4388, 0, 3, 701, 3590, 1016, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4472, 0, 3, 722, 3653, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4556, 0, 3, 743, 3716, 1072, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4640, 0, 3, 806, 3905, 1128, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4724, 0, 3, 827, 3968, 1156, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4808, 0, 3, 848, 4031, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4892, 0, 3, 869, 4094, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4976, 0, 3, 890, 4157, 1240, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5060, 0, 3, 932, 4220, 1268, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5168, 0, 3, 960, 4304, 1304, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5276, 0, 3, 988, 4388, 1340, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5384, 0, 3, 1016, 4472, 1376, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5492, 0, 3, 1044, 4556, 1412, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5600, 0, 3, 1100, 4640, 1448, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5708, 0, 3, 1128, 4724, 1484, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5816, 0, 3, 1156, 4808, 1520, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5924, 0, 3, 1184, 4892, 1556, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6032, 0, 3, 1212, 4976, 1592, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6140, 3, 9, 10, 1631, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6146, 3, 10, 11, 1634, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6152, 3, 11, 12, 1637, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6158, 3, 12, 13, 1640, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6164, 3, 13, 14, 1643, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6170, 3, 14, 15, 1646, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6176, 3, 15, 16, 1649, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6182, 3, 16, 17, 1652, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6188, 3, 22, 23, 1658, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6194, 3, 23, 24, 1661, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6200, 3, 24, 25, 1664, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6206, 3, 25, 26, 1667, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6212, 3, 26, 27, 1670, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6218, 3, 27, 28, 1673, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6224, 3, 28, 29, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6230, 3, 29, 30, 1679, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6236, 0, 3, 1628, 6140, 1691, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6254, 0, 3, 1631, 6146, 1700, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6272, 0, 3, 1634, 6152, 1709, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6290, 0, 3, 1637, 6158, 1718, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6308, 0, 3, 1640, 6164, 1727, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6326, 0, 3, 1643, 6170, 1736, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6344, 0, 3, 1646, 6176, 1745, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6362, 0, 3, 1649, 6182, 1754, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6380, 0, 3, 1655, 6188, 1772, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6398, 0, 3, 1658, 6194, 1781, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6416, 0, 3, 1661, 6200, 1790, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6434, 0, 3, 1664, 6206, 1799, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6452, 0, 3, 1667, 6212, 1808, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6470, 0, 3, 1670, 6218, 1817, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6488, 0, 3, 1673, 6224, 1826, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6506, 0, 3, 1676, 6230, 1835, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6524, 0, 3, 1682, 6236, 98, 104, 1862,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6560, 0, 3, 1691, 6254, 104, 110, 1880,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6596, 0, 3, 1700, 6272, 110, 116, 1898,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6632, 0, 3, 1709, 6290, 116, 122, 1916,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6668, 0, 3, 1718, 6308, 122, 128, 1934,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6704, 0, 3, 1727, 6326, 128, 134, 1952,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6740, 0, 3, 1736, 6344, 134, 140, 1970,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6776, 0, 3, 1745, 6362, 140, 146, 1988,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6812, 0, 3, 1763, 6380, 158, 164, 2024,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6848, 0, 3, 1772, 6398, 164, 170, 2042,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6884, 0, 3, 1781, 6416, 170, 176, 2060,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6920, 0, 3, 1790, 6434, 176, 182, 2078,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6956, 0, 3, 1799, 6452, 182, 188, 2096,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6992, 0, 3, 1808, 6470, 188, 194, 2114,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7028, 0, 3, 1817, 6488, 194, 200, 2132,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7064, 0, 3, 1826, 6506, 200, 206, 2150,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7100, 0, 3, 1862, 6560, 218, 228, 2228,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7160, 0, 3, 1880, 6596, 228, 238, 2258,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7220, 0, 3, 1898, 6632, 238, 248, 2288,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7280, 0, 3, 1916, 6668, 248, 258, 2318,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7340, 0, 3, 1934, 6704, 258, 268, 2348,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7400, 0, 3, 1952, 6740, 268, 278, 2378,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7460, 0, 3, 1970, 6776, 278, 288, 2408,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7520, 0, 3, 2024, 6848, 308, 318, 2498,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7580, 0, 3, 2042, 6884, 318, 328, 2528,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7640, 0, 3, 2060, 6920, 328, 338, 2558,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7700, 0, 3, 2078, 6956, 338, 348, 2588,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7760, 0, 3, 2096, 6992, 348, 358, 2618,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7820, 0, 3, 2114, 7028, 358, 368, 2648,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7880, 0, 3, 2132, 7064, 368, 378, 2678,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7940, 0, 3, 6524, 6560, 2228, 7160, 398,
                                                 413, 2753, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8030, 0, 3, 6560, 6596, 2258, 7220, 413,
                                                 428, 2798, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8120, 0, 3, 6596, 6632, 2288, 7280, 428,
                                                 443, 2843, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8210, 0, 3, 6632, 6668, 2318, 7340, 443,
                                                 458, 2888, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8300, 0, 3, 6668, 6704, 2348, 7400, 458,
                                                 473, 2933, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8390, 0, 3, 6704, 6740, 2378, 7460, 473,
                                                 488, 2978, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8480, 0, 3, 6812, 6848, 2498, 7580, 518,
                                                 533, 3068, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8570, 0, 3, 6848, 6884, 2528, 7640, 533,
                                                 548, 3113, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8660, 0, 3, 6884, 6920, 2558, 7700, 548,
                                                 563, 3158, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8750, 0, 3, 6920, 6956, 2588, 7760, 563,
                                                 578, 3203, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8840, 0, 3, 6956, 6992, 2618, 7820, 578,
                                                 593, 3248, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8930, 0, 3, 6992, 7028, 2648, 7880, 593,
                                                 608, 3293, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9020, 0, 3, 7100, 7160, 2753, 8030, 638,
                                                 659, 3464, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9146, 0, 3, 7160, 7220, 2798, 8120, 659,
                                                 680, 3527, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9272, 0, 3, 7220, 7280, 2843, 8210, 680,
                                                 701, 3590, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9398, 0, 3, 7280, 7340, 2888, 8300, 701,
                                                 722, 3653, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9524, 0, 3, 7340, 7400, 2933, 8390, 722,
                                                 743, 3716, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9650, 0, 3, 7520, 7580, 3068, 8570, 785,
                                                 806, 3905, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9776, 0, 3, 7580, 7640, 3113, 8660, 806,
                                                 827, 3968, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9902, 0, 3, 7640, 7700, 3158, 8750, 827,
                                                 848, 4031, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10028, 0, 3, 7700, 7760, 3203, 8840,
                                                 848, 869, 4094, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10154, 0, 3, 7760, 7820, 3248, 8930,
                                                 869, 890, 4157, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10280, 0, 3, 7940, 8030, 3464, 9146,
                                                 932, 960, 4304, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10448, 0, 3, 8030, 8120, 3527, 9272,
                                                 960, 988, 4388, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10616, 0, 3, 8120, 8210, 3590, 9398,
                                                 988, 1016, 4472, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10784, 0, 3, 8210, 8300, 3653, 9524,
                                                 1016, 1044, 4556, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10952, 0, 3, 8480, 8570, 3905, 9776,
                                                 1100, 1128, 4724, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11120, 0, 3, 8570, 8660, 3968, 9902,
                                                 1128, 1156, 4808, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11288, 0, 3, 8660, 8750, 4031, 10028,
                                                 1156, 1184, 4892, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11456, 0, 3, 8750, 8840, 4094, 10154,
                                                 1184, 1212, 4976, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11624, 0, 3, 9020, 9146, 4304, 10448,
                                                 1268, 1304, 5276, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11840, 0, 3, 9146, 9272, 4388, 10616,
                                                 1304, 1340, 5384, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12056, 0, 3, 9272, 9398, 4472, 10784,
                                                 1340, 1376, 5492, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12272, 0, 3, 9650, 9776, 4724, 11120,
                                                 1448, 1484, 5816, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12488, 0, 3, 9776, 9902, 4808, 11288,
                                                 1484, 1520, 5924, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12704, 0, 3, 9902, 10028, 4892, 11456,
                                                 1520, 1556, 6032, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12920, 3, 1628, 1631, 6146, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12930, 3, 1631, 1634, 6152, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12940, 3, 1634, 1637, 6158, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12950, 3, 1637, 1640, 6164, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12960, 3, 1640, 1643, 6170, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12970, 3, 1643, 1646, 6176, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12980, 3, 1646, 1649, 6182, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12990, 3, 1655, 1658, 6194, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13000, 3, 1658, 1661, 6200, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13010, 3, 1661, 1664, 6206, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13020, 3, 1664, 1667, 6212, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13030, 3, 1667, 1670, 6218, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13040, 3, 1670, 1673, 6224, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13050, 3, 1673, 1676, 6230, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 13060, 0, 3, 6140, 12920, 6254, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13090, 0, 3, 6146, 12930, 6272, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13120, 0, 3, 6152, 12940, 6290, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13150, 0, 3, 6158, 12950, 6308, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13180, 0, 3, 6164, 12960, 6326, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13210, 0, 3, 6170, 12970, 6344, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13240, 0, 3, 6176, 12980, 6362, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13270, 0, 3, 6188, 12990, 6398, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13300, 0, 3, 6194, 13000, 6416, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13330, 0, 3, 6200, 13010, 6434, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13360, 0, 3, 6206, 13020, 6452, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13390, 0, 3, 6212, 13030, 6470, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13420, 0, 3, 6218, 13040, 6488, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13450, 0, 3, 6224, 13050, 6506, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 13480, 0, 3, 6236, 13060, 1844, 1862,
                                                 6560, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13540, 0, 3, 6254, 13090, 1862, 1880,
                                                 6596, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13600, 0, 3, 6272, 13120, 1880, 1898,
                                                 6632, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13660, 0, 3, 6290, 13150, 1898, 1916,
                                                 6668, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13720, 0, 3, 6308, 13180, 1916, 1934,
                                                 6704, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13780, 0, 3, 6326, 13210, 1934, 1952,
                                                 6740, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13840, 0, 3, 6344, 13240, 1952, 1970,
                                                 6776, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13900, 0, 3, 6380, 13270, 2006, 2024,
                                                 6848, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13960, 0, 3, 6398, 13300, 2024, 2042,
                                                 6884, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14020, 0, 3, 6416, 13330, 2042, 2060,
                                                 6920, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14080, 0, 3, 6434, 13360, 2060, 2078,
                                                 6956, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14140, 0, 3, 6452, 13390, 2078, 2096,
                                                 6992, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14200, 0, 3, 6470, 13420, 2096, 2114,
                                                 7028, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14260, 0, 3, 6488, 13450, 2114, 2132,
                                                 7064, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14320, 0, 3, 6524, 13480, 2168, 2198,
                                                 7100, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14420, 0, 3, 6560, 13540, 2198, 2228,
                                                 7160, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14520, 0, 3, 6596, 13600, 2228, 2258,
                                                 7220, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14620, 0, 3, 6632, 13660, 2258, 2288,
                                                 7280, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14720, 0, 3, 6668, 13720, 2288, 2318,
                                                 7340, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14820, 0, 3, 6704, 13780, 2318, 2348,
                                                 7400, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14920, 0, 3, 6740, 13840, 2348, 2378,
                                                 7460, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15020, 0, 3, 6812, 13900, 2438, 2468,
                                                 7520, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15120, 0, 3, 6848, 13960, 2468, 2498,
                                                 7580, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15220, 0, 3, 6884, 14020, 2498, 2528,
                                                 7640, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15320, 0, 3, 6920, 14080, 2528, 2558,
                                                 7700, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15420, 0, 3, 6956, 14140, 2558, 2588,
                                                 7760, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15520, 0, 3, 6992, 14200, 2588, 2618,
                                                 7820, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15620, 0, 3, 7028, 14260, 2618, 2648,
                                                 7880, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15720, 0, 3, 13480, 13540, 7160, 14520,
                                                 2708, 2753, 8030, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15870, 0, 3, 13540, 13600, 7220, 14620,
                                                 2753, 2798, 8120, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16020, 0, 3, 13600, 13660, 7280, 14720,
                                                 2798, 2843, 8210, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16170, 0, 3, 13660, 13720, 7340, 14820,
                                                 2843, 2888, 8300, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16320, 0, 3, 13720, 13780, 7400, 14920,
                                                 2888, 2933, 8390, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16470, 0, 3, 13900, 13960, 7580, 15220,
                                                 3023, 3068, 8570, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16620, 0, 3, 13960, 14020, 7640, 15320,
                                                 3068, 3113, 8660, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16770, 0, 3, 14020, 14080, 7700, 15420,
                                                 3113, 3158, 8750, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16920, 0, 3, 14080, 14140, 7760, 15520,
                                                 3158, 3203, 8840, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17070, 0, 3, 14140, 14200, 7820, 15620,
                                                 3203, 3248, 8930, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17220, 0, 3, 14320, 14420, 7940, 15720,
                                                 3338, 3401, 9020, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17430, 0, 3, 14420, 14520, 8030, 15870,
                                                 3401, 3464, 9146, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17640, 0, 3, 14520, 14620, 8120, 16020,
                                                 3464, 3527, 9272, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17850, 0, 3, 14620, 14720, 8210, 16170,
                                                 3527, 3590, 9398, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18060, 0, 3, 14720, 14820, 8300, 16320,
                                                 3590, 3653, 9524, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18270, 0, 3, 15020, 15120, 8480, 16470,
                                                 3779, 3842, 9650, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18480, 0, 3, 15120, 15220, 8570, 16620,
                                                 3842, 3905, 9776, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18690, 0, 3, 15220, 15320, 8660, 16770,
                                                 3905, 3968, 9902, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18900, 0, 3, 15320, 15420, 8750, 16920,
                                                 3968, 4031, 10028, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19110, 0, 3, 15420, 15520, 8840, 17070,
                                                 4031, 4094, 10154, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19320, 0, 3, 15720, 15870, 9146, 17640,
                                                 4220, 4304, 10448, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19600, 0, 3, 15870, 16020, 9272, 17850,
                                                 4304, 4388, 10616, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19880, 0, 3, 16020, 16170, 9398, 18060,
                                                 4388, 4472, 10784, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20160, 0, 3, 16470, 16620, 9776, 18690,
                                                 4640, 4724, 11120, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20440, 0, 3, 16620, 16770, 9902, 18900,
                                                 4724, 4808, 11288, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20720, 0, 3, 16770, 16920, 10028, 19110,
                                                 4808, 4892, 11456, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 21000, 0, 3, 17220, 17430, 10280, 19320,
                                                 5060, 5168, 11624, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 21360, 0, 3, 17430, 17640, 10448, 19600,
                                                 5168, 5276, 11840, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 21720, 0, 3, 17640, 17850, 10616, 19880,
                                                 5276, 5384, 12056, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22080, 0, 3, 18270, 18480, 10952, 20160,
                                                 5600, 5708, 12272, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22440, 0, 3, 18480, 18690, 11120, 20440,
                                                 5708, 5816, 12488, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22800, 0, 3, 18690, 18900, 11288, 20720,
                                                 5816, 5924, 12704, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23160, 3, 6140, 6146, 12930, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23175, 3, 6146, 6152, 12940, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23190, 3, 6152, 6158, 12950, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23205, 3, 6158, 6164, 12960, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23220, 3, 6164, 6170, 12970, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23235, 3, 6170, 6176, 12980, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23250, 3, 6188, 6194, 13000, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23265, 3, 6194, 6200, 13010, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23280, 3, 6200, 6206, 13020, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23295, 3, 6206, 6212, 13030, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23310, 3, 6212, 6218, 13040, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23325, 3, 6218, 6224, 13050, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 23340, 0, 3, 12920, 23160, 13090, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23385, 0, 3, 12930, 23175, 13120, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23430, 0, 3, 12940, 23190, 13150, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23475, 0, 3, 12950, 23205, 13180, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23520, 0, 3, 12960, 23220, 13210, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23565, 0, 3, 12970, 23235, 13240, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23610, 0, 3, 12990, 23250, 13300, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23655, 0, 3, 13000, 23265, 13330, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23700, 0, 3, 13010, 23280, 13360, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23745, 0, 3, 13020, 23295, 13390, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23790, 0, 3, 13030, 23310, 13420, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23835, 0, 3, 13040, 23325, 13450, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 23880, 0, 3, 13060, 23340, 6524, 6560,
                                                 13540, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 23970, 0, 3, 13090, 23385, 6560, 6596,
                                                 13600, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24060, 0, 3, 13120, 23430, 6596, 6632,
                                                 13660, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24150, 0, 3, 13150, 23475, 6632, 6668,
                                                 13720, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24240, 0, 3, 13180, 23520, 6668, 6704,
                                                 13780, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24330, 0, 3, 13210, 23565, 6704, 6740,
                                                 13840, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24420, 0, 3, 13270, 23610, 6812, 6848,
                                                 13960, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24510, 0, 3, 13300, 23655, 6848, 6884,
                                                 14020, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24600, 0, 3, 13330, 23700, 6884, 6920,
                                                 14080, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24690, 0, 3, 13360, 23745, 6920, 6956,
                                                 14140, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24780, 0, 3, 13390, 23790, 6956, 6992,
                                                 14200, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24870, 0, 3, 13420, 23835, 6992, 7028,
                                                 14260, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24960, 0, 3, 13540, 23970, 7100, 7160,
                                                 14520, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25110, 0, 3, 13600, 24060, 7160, 7220,
                                                 14620, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25260, 0, 3, 13660, 24150, 7220, 7280,
                                                 14720, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25410, 0, 3, 13720, 24240, 7280, 7340,
                                                 14820, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25560, 0, 3, 13780, 24330, 7340, 7400,
                                                 14920, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25710, 0, 3, 13960, 24510, 7520, 7580,
                                                 15220, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25860, 0, 3, 14020, 24600, 7580, 7640,
                                                 15320, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 26010, 0, 3, 14080, 24690, 7640, 7700,
                                                 15420, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 26160, 0, 3, 14140, 24780, 7700, 7760,
                                                 15520, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 26310, 0, 3, 14200, 24870, 7760, 7820,
                                                 15620, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26460, 0, 3, 23880, 23970, 14520, 25110,
                                                 7940, 8030, 15870, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26685, 0, 3, 23970, 24060, 14620, 25260,
                                                 8030, 8120, 16020, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26910, 0, 3, 24060, 24150, 14720, 25410,
                                                 8120, 8210, 16170, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27135, 0, 3, 24150, 24240, 14820, 25560,
                                                 8210, 8300, 16320, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27360, 0, 3, 24420, 24510, 15220, 25860,
                                                 8480, 8570, 16620, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27585, 0, 3, 24510, 24600, 15320, 26010,
                                                 8570, 8660, 16770, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27810, 0, 3, 24600, 24690, 15420, 26160,
                                                 8660, 8750, 16920, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 28035, 0, 3, 24690, 24780, 15520, 26310,
                                                 8750, 8840, 17070, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28260, 0, 3, 24960, 25110, 15870, 26685,
                                                 9020, 9146, 17640, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28575, 0, 3, 25110, 25260, 16020, 26910,
                                                 9146, 9272, 17850, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28890, 0, 3, 25260, 25410, 16170, 27135,
                                                 9272, 9398, 18060, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29205, 0, 3, 25710, 25860, 16620, 27585,
                                                 9650, 9776, 18690, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29520, 0, 3, 25860, 26010, 16770, 27810,
                                                 9776, 9902, 18900, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29835, 0, 3, 26010, 26160, 16920, 28035,
                                                 9902, 10028, 19110, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 30150, 0, 3, 26460, 26685, 17640, 28575,
                                                 10280, 10448, 19600, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 30570, 0, 3, 26685, 26910, 17850, 28890,
                                                 10448, 10616, 19880, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 30990, 0, 3, 27360, 27585, 18690, 29520,
                                                 10952, 11120, 20440, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 31410, 0, 3, 27585, 27810, 18900, 29835,
                                                 11120, 11288, 20720, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 31830, 0, 3, 28260, 28575, 19600, 30570,
                                                 11624, 11840, 21720, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 32370, 0, 3, 29205, 29520, 20440, 31410,
                                                 12272, 12488, 22800, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 32910, 3, 12920, 12930, 23175, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 32931, 3, 12930, 12940, 23190, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 32952, 3, 12940, 12950, 23205, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 32973, 3, 12950, 12960, 23220, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 32994, 3, 12960, 12970, 23235, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33015, 3, 12990, 13000, 23265, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33036, 3, 13000, 13010, 23280, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33057, 3, 13010, 13020, 23295, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33078, 3, 13020, 13030, 23310, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 33099, 3, 13030, 13040, 23325, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 33120, 0, 3, 23160, 32910, 23385, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33183, 0, 3, 23175, 32931, 23430, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33246, 0, 3, 23190, 32952, 23475, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33309, 0, 3, 23205, 32973, 23520, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33372, 0, 3, 23220, 32994, 23565, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33435, 0, 3, 23250, 33015, 23655, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33498, 0, 3, 23265, 33036, 23700, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33561, 0, 3, 23280, 33057, 23745, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33624, 0, 3, 23295, 33078, 23790, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 33687, 0, 3, 23310, 33099, 23835, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 33750, 0, 3, 23340, 33120, 13480, 13540,
                                                 23970, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 33876, 0, 3, 23385, 33183, 13540, 13600,
                                                 24060, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34002, 0, 3, 23430, 33246, 13600, 13660,
                                                 24150, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34128, 0, 3, 23475, 33309, 13660, 13720,
                                                 24240, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34254, 0, 3, 23520, 33372, 13720, 13780,
                                                 24330, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34380, 0, 3, 23610, 33435, 13900, 13960,
                                                 24510, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34506, 0, 3, 23655, 33498, 13960, 14020,
                                                 24600, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34632, 0, 3, 23700, 33561, 14020, 14080,
                                                 24690, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34758, 0, 3, 23745, 33624, 14080, 14140,
                                                 24780, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34884, 0, 3, 23790, 33687, 14140, 14200,
                                                 24870, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35010, 0, 3, 23880, 33750, 14320, 14420,
                                                 24960, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35220, 0, 3, 23970, 33876, 14420, 14520,
                                                 25110, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35430, 0, 3, 24060, 34002, 14520, 14620,
                                                 25260, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35640, 0, 3, 24150, 34128, 14620, 14720,
                                                 25410, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35850, 0, 3, 24240, 34254, 14720, 14820,
                                                 25560, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36060, 0, 3, 24420, 34380, 15020, 15120,
                                                 25710, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36270, 0, 3, 24510, 34506, 15120, 15220,
                                                 25860, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36480, 0, 3, 24600, 34632, 15220, 15320,
                                                 26010, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36690, 0, 3, 24690, 34758, 15320, 15420,
                                                 26160, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36900, 0, 3, 24780, 34884, 15420, 15520,
                                                 26310, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37110, 0, 3, 33750, 33876, 25110, 35430,
                                                 15720, 15870, 26685, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37425, 0, 3, 33876, 34002, 25260, 35640,
                                                 15870, 16020, 26910, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37740, 0, 3, 34002, 34128, 25410, 35850,
                                                 16020, 16170, 27135, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38055, 0, 3, 34380, 34506, 25860, 36480,
                                                 16470, 16620, 27585, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38370, 0, 3, 34506, 34632, 26010, 36690,
                                                 16620, 16770, 27810, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38685, 0, 3, 34632, 34758, 26160, 36900,
                                                 16770, 16920, 28035, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39000, 0, 3, 35010, 35220, 26460, 37110,
                                                 17220, 17430, 28260, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39441, 0, 3, 35220, 35430, 26685, 37425,
                                                 17430, 17640, 28575, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39882, 0, 3, 35430, 35640, 26910, 37740,
                                                 17640, 17850, 28890, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 40323, 0, 3, 36060, 36270, 27360, 38055,
                                                 18270, 18480, 29205, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 40764, 0, 3, 36270, 36480, 27585, 38370,
                                                 18480, 18690, 29520, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 41205, 0, 3, 36480, 36690, 27810, 38685,
                                                 18690, 18900, 29835, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 41646, 0, 3, 37110, 37425, 28575, 39882,
                                                 19320, 19600, 30570, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 42234, 0, 3, 38055, 38370, 29520, 41205,
                                                 20160, 20440, 31410, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 42822, 0, 3, 39000, 39441, 30150, 41646,
                                                 21000, 21360, 31830, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 43578, 0, 3, 40323, 40764, 30990, 42234,
                                                 22080, 22440, 32370, ncols, alpha, beta, p);

            simdgeo::geom_i_x(buffer, 44334, 40323, 43578, 1, 21, ncols, alpha);

            simdgeo::geom_i_y(buffer, 44922, 40323, 43578, 1, 21, ncols, alpha);

            simdgeo::geom_i_z(buffer, 45510, 40323, 43578, 1, 21, ncols, alpha);

            simdgeo::geom_i_x(buffer, 46098, 39000, 42822, 1, 21, ncols, alpha);

            simdgeo::geom_i_y(buffer, 46686, 39000, 42822, 1, 21, ncols, alpha);

            simdgeo::geom_i_z(buffer, 47274, 39000, 42822, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 47862, 46098, 1764, ncols);

            simdfunc::contract_primitives(buffer, 49626, 44334, 1764, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 51390, 49626, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 51390, 11, nmax);

    simdtrf::transform_h_inner(buffer, 51390, 50214, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 143 * nvalues, nvalues, buffer, 51390, 11, nmax);

    simdtrf::transform_h_inner(buffer, 51390, 50802, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 286 * nvalues, nvalues, buffer, 51390, 11, nmax);

    simdtrf::transform_h_inner(buffer, 51390, 47862, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 429 * nvalues, nvalues, buffer, 51390, 11, nmax);

    simdtrf::transform_h_inner(buffer, 51390, 48450, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 572 * nvalues, nvalues, buffer, 51390, 11, nmax);

    simdtrf::transform_h_inner(buffer, 51390, 49038, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 715 * nvalues, nvalues, buffer, 51390, 11, nmax);
}

}  // namespace simdt2ceri
