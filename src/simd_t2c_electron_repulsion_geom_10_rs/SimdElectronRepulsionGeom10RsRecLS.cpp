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


#include "SimdElectronRepulsionGeom10RsRecLS.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ls_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ls_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1785, 1515, 270, nvalues);

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

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_erf_boys_function(buffer, coordinates, 3, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 13, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 4, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 5, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 6, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 22, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 4, 5, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 5, 6, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 6, 7, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 7, 8, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 8, 9, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 9, 10, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 10, 11, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 119, 0, 14, 15, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 125, 0, 15, 16, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 131, 0, 16, 17, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 137, 0, 17, 18, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 143, 0, 18, 19, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 149, 0, 19, 20, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 155, 0, 20, 21, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 161, 0, 23, 26, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 171, 0, 26, 29, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 181, 0, 29, 32, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 191, 0, 32, 35, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 201, 0, 35, 38, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 211, 0, 38, 41, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 221, 0, 41, 44, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 231, 0, 50, 53, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 241, 0, 53, 56, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 251, 0, 56, 59, 131, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 261, 0, 59, 62, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 271, 0, 62, 65, 143, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 281, 0, 65, 68, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 291, 0, 68, 71, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 301, 0, 77, 83, 181, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 316, 0, 83, 89, 191, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 331, 0, 89, 95, 201, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 346, 0, 95, 101, 211, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 361, 0, 101, 107, 221, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 376, 0, 119, 125, 251, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 391, 0, 125, 131, 261, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 406, 0, 131, 137, 271, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 421, 0, 137, 143, 281, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 436, 0, 143, 149, 291, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 451, 0, 161, 171, 301, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 472, 0, 171, 181, 316, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 493, 0, 181, 191, 331, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 514, 0, 191, 201, 346, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 535, 0, 201, 211, 361, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 556, 0, 231, 241, 376, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 577, 0, 241, 251, 391, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 598, 0, 251, 261, 406, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 619, 0, 261, 271, 421, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 640, 0, 271, 281, 436, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 661, 0, 301, 316, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 689, 0, 316, 331, 514, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 717, 0, 331, 346, 535, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 745, 0, 376, 391, 598, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 773, 0, 391, 406, 619, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 801, 0, 406, 421, 640, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 829, 0, 451, 472, 661, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 865, 0, 472, 493, 689, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 901, 0, 493, 514, 717, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 937, 0, 556, 577, 745, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 973, 0, 577, 598, 773, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1009, 0, 598, 619, 801, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1045, 0, 661, 689, 901, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1090, 0, 745, 773, 1009, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1135, 0, 829, 865, 1045, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1190, 0, 937, 973, 1090, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_l_x(buffer, 1245, 937, 1190, 1, 1, ncols, alpha);

            simdgeo::geom_l_y(buffer, 1290, 937, 1190, 1, 1, ncols, alpha);

            simdgeo::geom_l_z(buffer, 1335, 937, 1190, 1, 1, ncols, alpha);

            simdgeo::geom_l_x(buffer, 1380, 829, 1135, 1, 1, ncols, alpha);

            simdgeo::geom_l_y(buffer, 1425, 829, 1135, 1, 1, ncols, alpha);

            simdgeo::geom_l_z(buffer, 1470, 829, 1135, 1, 1, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1515, 1380, 135, ncols);

            simdfunc::contract_primitives(buffer, 1650, 1245, 135, ncols);
        }
    }

    simdtrf::transform_l_outer(values, nvalues, buffer, 1650, 1, nmax);

    simdtrf::transform_l_outer(values + 17 * nvalues, nvalues, buffer, 1695, 1, nmax);

    simdtrf::transform_l_outer(values + 34 * nvalues, nvalues, buffer, 1740, 1, nmax);

    simdtrf::transform_l_outer(values + 51 * nvalues, nvalues, buffer, 1515, 1, nmax);

    simdtrf::transform_l_outer(values + 68 * nvalues, nvalues, buffer, 1560, 1, nmax);

    simdtrf::transform_l_outer(values + 85 * nvalues, nvalues, buffer, 1605, 1, nmax);
}

}  // namespace simdt2ceri
