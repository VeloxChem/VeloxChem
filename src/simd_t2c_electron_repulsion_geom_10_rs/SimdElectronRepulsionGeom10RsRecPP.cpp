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


#include "SimdElectronRepulsionGeom10RsRecPP.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 215, 152, 54, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3}, ncols, fj, mu,
                                                omega);

            simdfunc::compute_boys_function(buffer, coordinates, 10, {1, 2, 3}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 26, 0, 7, 8, 17, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 32, 0, 11, 12, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 38, 3, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 41, 3, 11, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 44, 3, 8, 17, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 53, 3, 12, 23, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 62, 0, 3, 14, 44, 26, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 80, 0, 3, 20, 53, 32, ncols, p);

            simdgeo::geom_p_x(buffer, 98, 41, 80, 1, 3, ncols, alpha);

            simdgeo::geom_p_y(buffer, 107, 41, 80, 1, 3, ncols, alpha);

            simdgeo::geom_p_z(buffer, 116, 41, 80, 1, 3, ncols, alpha);

            simdgeo::geom_p_x(buffer, 125, 38, 62, 1, 3, ncols, alpha);

            simdgeo::geom_p_y(buffer, 134, 38, 62, 1, 3, ncols, alpha);

            simdgeo::geom_p_z(buffer, 143, 38, 62, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 152, 125, 27, ncols);

            simdfunc::contract_primitives(buffer, 179, 98, 27, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 206, 179, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 206, 3, nmax);

    simdtrf::transform_p_inner(buffer, 206, 188, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 9 * nvalues, nvalues, buffer, 206, 3, nmax);

    simdtrf::transform_p_inner(buffer, 206, 197, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 18 * nvalues, nvalues, buffer, 206, 3, nmax);

    simdtrf::transform_p_inner(buffer, 206, 152, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 27 * nvalues, nvalues, buffer, 206, 3, nmax);

    simdtrf::transform_p_inner(buffer, 206, 161, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 36 * nvalues, nvalues, buffer, 206, 3, nmax);

    simdtrf::transform_p_inner(buffer, 206, 170, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 45 * nvalues, nvalues, buffer, 206, 3, nmax);
}

}  // namespace simdt2ceri
