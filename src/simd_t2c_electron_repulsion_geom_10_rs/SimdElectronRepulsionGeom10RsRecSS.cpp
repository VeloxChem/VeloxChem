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


#include "SimdElectronRepulsionGeom10RsRecSS.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryS1.hpp"
#include "SimdTransformS.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ss_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ss_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 25, 19, 6, nvalues);

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

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_erf_boys_function(buffer, coordinates, 3, {1}, ncols, fj, mu,
                                                omega);

            simdfunc::compute_boys_function(buffer, coordinates, 5, {1}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 7, 0, 4, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 10, 0, 6, ncols);

            simdgeo::geom_s_x(buffer, 13, 10, 1, 1, ncols, alpha);

            simdgeo::geom_s_y(buffer, 14, 10, 1, 1, ncols, alpha);

            simdgeo::geom_s_z(buffer, 15, 10, 1, 1, ncols, alpha);

            simdgeo::geom_s_x(buffer, 16, 7, 1, 1, ncols, alpha);

            simdgeo::geom_s_y(buffer, 17, 7, 1, 1, ncols, alpha);

            simdgeo::geom_s_z(buffer, 18, 7, 1, 1, ncols, alpha);

            simdfunc::contract_primitives(buffer, 19, 16, 3, ncols);

            simdfunc::contract_primitives(buffer, 22, 13, 3, ncols);
        }
    }

    simdtrf::transform_s_outer(values, nvalues, buffer, 22, 1, nmax);

    simdtrf::transform_s_outer(values + 1 * nvalues, nvalues, buffer, 23, 1, nmax);

    simdtrf::transform_s_outer(values + 2 * nvalues, nvalues, buffer, 24, 1, nmax);

    simdtrf::transform_s_outer(values + 3 * nvalues, nvalues, buffer, 19, 1, nmax);

    simdtrf::transform_s_outer(values + 4 * nvalues, nvalues, buffer, 20, 1, nmax);

    simdtrf::transform_s_outer(values + 5 * nvalues, nvalues, buffer, 21, 1, nmax);
}

}  // namespace simdt2ceri
