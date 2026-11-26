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


#include "SimdElectronRepulsionGeom10RsRecDP.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_dp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 444, 318, 108, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4}, ncols, fj,
                                                mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 11, {1, 2, 3, 4}, ncols, fj,
                                            mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 34, 0, 7, 8, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 8, 9, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 12, 13, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 13, 14, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 58, 0, 16, 19, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 68, 0, 25, 28, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_pp_electron_repulsion_0(buffer, 78, 3, 7, 16, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 87, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 96, 3, 12, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 105, 3, 14, 31, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 114, 0, 3, 19, 87, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 132, 0, 3, 28, 105, 52, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 150, 0, 3, 34, 114, 58, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 180, 0, 3, 46, 132, 68, ncols, p);

            simdgeo::geom_d_x(buffer, 210, 96, 180, 1, 3, ncols, alpha);

            simdgeo::geom_d_y(buffer, 228, 96, 180, 1, 3, ncols, alpha);

            simdgeo::geom_d_z(buffer, 246, 96, 180, 1, 3, ncols, alpha);

            simdgeo::geom_d_x(buffer, 264, 78, 150, 1, 3, ncols, alpha);

            simdgeo::geom_d_y(buffer, 282, 78, 150, 1, 3, ncols, alpha);

            simdgeo::geom_d_z(buffer, 300, 78, 150, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 318, 264, 54, ncols);

            simdfunc::contract_primitives(buffer, 372, 210, 54, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 426, 372, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 426, 3, nmax);

    simdtrf::transform_p_inner(buffer, 426, 390, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 15 * nvalues, nvalues, buffer, 426, 3, nmax);

    simdtrf::transform_p_inner(buffer, 426, 408, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 30 * nvalues, nvalues, buffer, 426, 3, nmax);

    simdtrf::transform_p_inner(buffer, 426, 318, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 45 * nvalues, nvalues, buffer, 426, 3, nmax);

    simdtrf::transform_p_inner(buffer, 426, 336, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 60 * nvalues, nvalues, buffer, 426, 3, nmax);

    simdtrf::transform_p_inner(buffer, 426, 354, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 75 * nvalues, nvalues, buffer, 426, 3, nmax);
}

}  // namespace simdt2ceri
