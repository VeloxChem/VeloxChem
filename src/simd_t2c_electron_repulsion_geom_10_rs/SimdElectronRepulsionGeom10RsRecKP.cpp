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


#include "SimdElectronRepulsionGeom10RsRecKP.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_kp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_kp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3884, 3128, 648, nvalues);

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
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 16, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 7, 8, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 8, 9, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 9, 10, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 10, 11, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 11, 12, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 12, 13, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 13, 14, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 18, 19, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 19, 20, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 20, 21, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 21, 22, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 22, 23, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 23, 24, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 50, 53, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 53, 56, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 56, 59, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 59, 62, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 62, 65, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 65, 68, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 278, 0, 74, 80, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 293, 0, 80, 86, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 308, 0, 86, 92, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 323, 0, 92, 98, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 338, 0, 98, 104, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 353, 0, 116, 122, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 368, 0, 122, 128, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 383, 0, 128, 134, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 398, 0, 134, 140, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 413, 0, 140, 146, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 428, 0, 158, 168, 293, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 449, 0, 168, 178, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 470, 0, 178, 188, 323, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 491, 0, 188, 198, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 512, 0, 218, 228, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 533, 0, 228, 238, 383, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 554, 0, 238, 248, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 575, 0, 248, 258, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 596, 0, 278, 293, 449, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 624, 0, 293, 308, 470, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 652, 0, 308, 323, 491, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 680, 0, 353, 368, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 708, 0, 368, 383, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 736, 0, 383, 398, 575, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 764, 0, 428, 449, 624, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 800, 0, 449, 470, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 836, 0, 512, 533, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 872, 0, 533, 554, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 908, 0, 596, 624, 800, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 953, 0, 680, 708, 872, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 998, 3, 12, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1007, 3, 14, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1016, 3, 22, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1025, 3, 24, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1034, 0, 3, 38, 998, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1052, 0, 3, 44, 1007, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1070, 0, 3, 62, 1016, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1088, 0, 3, 68, 1025, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1106, 0, 3, 92, 1034, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1136, 0, 3, 104, 1052, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1166, 0, 3, 134, 1070, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1196, 0, 3, 146, 1088, 268, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1226, 0, 3, 178, 1106, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1271, 0, 3, 198, 1136, 338, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1316, 0, 3, 238, 1166, 383, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1361, 0, 3, 258, 1196, 413, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1406, 0, 3, 293, 1226, 449, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1469, 0, 3, 323, 1271, 491, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1532, 0, 3, 368, 1316, 533, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1595, 0, 3, 398, 1361, 575, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1658, 0, 3, 428, 1406, 596, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1742, 0, 3, 470, 1469, 652, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1826, 0, 3, 512, 1532, 680, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1910, 0, 3, 554, 1595, 736, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1994, 0, 3, 624, 1742, 800, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2102, 0, 3, 708, 1910, 872, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2210, 0, 3, 764, 1994, 908, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2345, 0, 3, 836, 2102, 953, ncols, p);

            simdgeo::geom_k_x(buffer, 2480, 1826, 2345, 1, 3, ncols, alpha);

            simdgeo::geom_k_y(buffer, 2588, 1826, 2345, 1, 3, ncols, alpha);

            simdgeo::geom_k_z(buffer, 2696, 1826, 2345, 1, 3, ncols, alpha);

            simdgeo::geom_k_x(buffer, 2804, 1658, 2210, 1, 3, ncols, alpha);

            simdgeo::geom_k_y(buffer, 2912, 1658, 2210, 1, 3, ncols, alpha);

            simdgeo::geom_k_z(buffer, 3020, 1658, 2210, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 3128, 2804, 324, ncols);

            simdfunc::contract_primitives(buffer, 3452, 2480, 324, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 3776, 3452, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 3776, 3, nmax);

    simdtrf::transform_p_inner(buffer, 3776, 3560, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 45 * nvalues, nvalues, buffer, 3776, 3, nmax);

    simdtrf::transform_p_inner(buffer, 3776, 3668, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 90 * nvalues, nvalues, buffer, 3776, 3, nmax);

    simdtrf::transform_p_inner(buffer, 3776, 3128, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 135 * nvalues, nvalues, buffer, 3776, 3, nmax);

    simdtrf::transform_p_inner(buffer, 3776, 3236, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 180 * nvalues, nvalues, buffer, 3776, 3, nmax);

    simdtrf::transform_p_inner(buffer, 3776, 3344, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 225 * nvalues, nvalues, buffer, 3776, 3, nmax);
}

}  // namespace simdt2ceri
