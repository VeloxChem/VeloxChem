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


#include "SimdElectronRepulsionGeom10RsRecGS.hpp"

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
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_gs_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_gs_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 393, 303, 90, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 3, {1, 2, 3, 4, 5}, ncols,
                                                fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 9, {1, 2, 3, 4, 5}, ncols, fj,
                                            mu);

            compute_prim_ps_electron_repulsion_0(buffer, 15, 0, 4, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 5, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 6, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 14, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 45, 0, 4, 5, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 5, 6, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 6, 7, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 10, 11, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 11, 12, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 12, 13, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 81, 0, 15, 18, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 91, 0, 18, 21, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 101, 0, 21, 24, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 111, 0, 30, 33, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 121, 0, 33, 36, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 131, 0, 36, 39, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 141, 0, 45, 51, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 156, 0, 63, 69, 131, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 171, 0, 81, 91, 141, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 192, 0, 111, 121, 156, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_g_x(buffer, 213, 111, 192, 1, 1, ncols, alpha);

            simdgeo::geom_g_y(buffer, 228, 111, 192, 1, 1, ncols, alpha);

            simdgeo::geom_g_z(buffer, 243, 111, 192, 1, 1, ncols, alpha);

            simdgeo::geom_g_x(buffer, 258, 81, 171, 1, 1, ncols, alpha);

            simdgeo::geom_g_y(buffer, 273, 81, 171, 1, 1, ncols, alpha);

            simdgeo::geom_g_z(buffer, 288, 81, 171, 1, 1, ncols, alpha);

            simdfunc::contract_primitives(buffer, 303, 258, 45, ncols);

            simdfunc::contract_primitives(buffer, 348, 213, 45, ncols);
        }
    }

    simdtrf::transform_g_outer(values, nvalues, buffer, 348, 1, nmax);

    simdtrf::transform_g_outer(values + 9 * nvalues, nvalues, buffer, 363, 1, nmax);

    simdtrf::transform_g_outer(values + 18 * nvalues, nvalues, buffer, 378, 1, nmax);

    simdtrf::transform_g_outer(values + 27 * nvalues, nvalues, buffer, 303, 1, nmax);

    simdtrf::transform_g_outer(values + 36 * nvalues, nvalues, buffer, 318, 1, nmax);

    simdtrf::transform_g_outer(values + 45 * nvalues, nvalues, buffer, 333, 1, nmax);
}

}  // namespace simdt2ceri
