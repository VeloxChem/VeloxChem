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


#include "SimdElectronRepulsionGeom10RsRecLP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_lp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_lp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5221, 4276, 810, nvalues);

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
                                                9, 10}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 17, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 27, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 7, 8, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 8, 9, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 9, 10, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 10, 11, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 11, 12, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 12, 13, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 13, 14, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 14, 15, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 20, 21, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 21, 22, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 22, 23, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 23, 24, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 24, 25, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 25, 26, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 28, 31, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 31, 34, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 34, 37, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 37, 40, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 40, 43, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 43, 46, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 46, 49, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 55, 58, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 58, 61, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 61, 64, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 64, 67, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 67, 70, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 70, 73, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 73, 76, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 318, 0, 82, 88, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 333, 0, 88, 94, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 348, 0, 94, 100, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 363, 0, 100, 106, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 378, 0, 106, 112, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 393, 0, 112, 118, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 408, 0, 130, 136, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 423, 0, 136, 142, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 438, 0, 142, 148, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 453, 0, 148, 154, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 468, 0, 154, 160, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 483, 0, 160, 166, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 498, 0, 178, 188, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 519, 0, 188, 198, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 540, 0, 198, 208, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 561, 0, 208, 218, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 582, 0, 218, 228, 393, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 603, 0, 248, 258, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 624, 0, 258, 268, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 645, 0, 268, 278, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 666, 0, 278, 288, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 687, 0, 288, 298, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 708, 0, 318, 333, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 736, 0, 333, 348, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 764, 0, 348, 363, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 792, 0, 363, 378, 582, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 820, 0, 408, 423, 624, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 848, 0, 423, 438, 645, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 876, 0, 438, 453, 666, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 904, 0, 453, 468, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 932, 0, 498, 519, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 968, 0, 519, 540, 764, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1004, 0, 540, 561, 792, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1040, 0, 603, 624, 848, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1076, 0, 624, 645, 876, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1112, 0, 645, 666, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1148, 0, 708, 736, 968, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1193, 0, 736, 764, 1004, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1238, 0, 820, 848, 1076, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1283, 0, 848, 876, 1112, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1328, 0, 932, 968, 1193, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1383, 0, 1040, 1076, 1283, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1438, 3, 13, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1447, 3, 15, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1456, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1465, 3, 26, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1474, 0, 3, 43, 1438, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1492, 0, 3, 49, 1447, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1510, 0, 3, 70, 1456, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1528, 0, 3, 76, 1465, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1546, 0, 3, 106, 1474, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1576, 0, 3, 118, 1492, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1606, 0, 3, 154, 1510, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1636, 0, 3, 166, 1528, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1666, 0, 3, 208, 1546, 363, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1711, 0, 3, 228, 1576, 393, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1756, 0, 3, 278, 1606, 453, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1801, 0, 3, 298, 1636, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1846, 0, 3, 348, 1666, 540, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1909, 0, 3, 378, 1711, 582, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1972, 0, 3, 438, 1756, 645, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2035, 0, 3, 468, 1801, 687, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2098, 0, 3, 519, 1846, 736, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2182, 0, 3, 561, 1909, 792, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2266, 0, 3, 624, 1972, 848, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2350, 0, 3, 666, 2035, 904, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2434, 0, 3, 708, 2098, 932, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2542, 0, 3, 764, 2182, 1004, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2650, 0, 3, 820, 2266, 1040, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2758, 0, 3, 876, 2350, 1112, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2866, 0, 3, 968, 2542, 1193, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3001, 0, 3, 1076, 2758, 1283, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 3136, 0, 3, 1148, 2866, 1328, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 3301, 0, 3, 1238, 3001, 1383, ncols,
                                                 p);

            simdgeo::geom_l_x(buffer, 3466, 2650, 3301, 1, 3, ncols, alpha);

            simdgeo::geom_l_y(buffer, 3601, 2650, 3301, 1, 3, ncols, alpha);

            simdgeo::geom_l_z(buffer, 3736, 2650, 3301, 1, 3, ncols, alpha);

            simdgeo::geom_l_x(buffer, 3871, 2434, 3136, 1, 3, ncols, alpha);

            simdgeo::geom_l_y(buffer, 4006, 2434, 3136, 1, 3, ncols, alpha);

            simdgeo::geom_l_z(buffer, 4141, 2434, 3136, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 4276, 3871, 405, ncols);

            simdfunc::contract_primitives(buffer, 4681, 3466, 405, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 5086, 4681, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 5086, 3, nmax);

    simdtrf::transform_p_inner(buffer, 5086, 4816, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 51 * nvalues, nvalues, buffer, 5086, 3, nmax);

    simdtrf::transform_p_inner(buffer, 5086, 4951, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 102 * nvalues, nvalues, buffer, 5086, 3, nmax);

    simdtrf::transform_p_inner(buffer, 5086, 4276, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 153 * nvalues, nvalues, buffer, 5086, 3, nmax);

    simdtrf::transform_p_inner(buffer, 5086, 4411, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 204 * nvalues, nvalues, buffer, 5086, 3, nmax);

    simdtrf::transform_p_inner(buffer, 5086, 4546, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 255 * nvalues, nvalues, buffer, 5086, 3, nmax);
}

}  // namespace simdt2ceri
