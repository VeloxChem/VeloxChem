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


#include "SimdElectronRepulsionGeom10RsRecKS.hpp"

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
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ks_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ks_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1287, 1071, 216, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 3, 8, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 13, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 6, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 22, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 4, 5, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 5, 6, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 6, 7, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 7, 8, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 8, 9, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 9, 10, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 10, 11, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 14, 15, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 15, 16, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 119, 0, 16, 17, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 125, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 131, 0, 18, 19, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 137, 0, 19, 20, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 143, 0, 20, 21, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 23, 26, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 26, 29, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 29, 32, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 32, 35, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 35, 38, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 44, 47, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 209, 0, 47, 50, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 219, 0, 50, 53, 131, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 229, 0, 53, 56, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 239, 0, 56, 59, 143, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 249, 0, 65, 71, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 264, 0, 71, 77, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 279, 0, 77, 83, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 294, 0, 83, 89, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 309, 0, 89, 95, 189, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 324, 0, 107, 113, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 339, 0, 113, 119, 209, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 354, 0, 119, 125, 219, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 369, 0, 125, 131, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 384, 0, 131, 137, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 399, 0, 149, 159, 279, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 420, 0, 159, 169, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 441, 0, 169, 179, 309, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 462, 0, 199, 209, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 483, 0, 209, 219, 369, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 504, 0, 219, 229, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 249, 264, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 553, 0, 264, 279, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 279, 294, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 324, 339, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 339, 354, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 354, 369, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 693, 0, 399, 420, 581, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 729, 0, 462, 483, 665, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 765, 0, 525, 553, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 810, 0, 609, 637, 729, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_k_x(buffer, 855, 609, 810, 1, 1, ncols, alpha);

            simdgeo::geom_k_y(buffer, 891, 609, 810, 1, 1, ncols, alpha);

            simdgeo::geom_k_z(buffer, 927, 609, 810, 1, 1, ncols, alpha);

            simdgeo::geom_k_x(buffer, 963, 525, 765, 1, 1, ncols, alpha);

            simdgeo::geom_k_y(buffer, 999, 525, 765, 1, 1, ncols, alpha);

            simdgeo::geom_k_z(buffer, 1035, 525, 765, 1, 1, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1071, 963, 108, ncols);

            simdfunc::contract_primitives(buffer, 1179, 855, 108, ncols);
        }
    }

    simdtrf::transform_k_outer(values, nvalues, buffer, 1179, 1, nmax);

    simdtrf::transform_k_outer(values + 15 * nvalues, nvalues, buffer, 1215, 1, nmax);

    simdtrf::transform_k_outer(values + 30 * nvalues, nvalues, buffer, 1251, 1, nmax);

    simdtrf::transform_k_outer(values + 45 * nvalues, nvalues, buffer, 1071, 1, nmax);

    simdtrf::transform_k_outer(values + 60 * nvalues, nvalues, buffer, 1107, 1, nmax);

    simdtrf::transform_k_outer(values + 75 * nvalues, nvalues, buffer, 1143, 1, nmax);
}

}  // namespace simdt2ceri
