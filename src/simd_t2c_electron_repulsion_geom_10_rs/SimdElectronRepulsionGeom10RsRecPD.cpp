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


#include "SimdElectronRepulsionGeom10RsRecPD.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pd_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 525, 402, 108, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 4, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 12, 4, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 36, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 42, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 13, 14, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 14, 15, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 15, 16, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 72, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 75, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 78, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 81, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 84, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 87, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 90, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 99, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 108, 3, 15, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 117, 3, 16, 33, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 126, 0, 3, 21, 99, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 144, 0, 3, 30, 117, 66, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 162, 3, 7, 8, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 168, 3, 9, 10, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 174, 3, 13, 14, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 180, 3, 15, 16, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 186, 0, 3, 75, 168, 99, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 204, 0, 3, 84, 180, 117, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 222, 0, 3, 90, 186, 36, 42, 126, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 258, 0, 3, 108, 204, 54, 60, 144, ncols,
                                                 alpha, beta, p);

            simdgeo::geom_p_x(buffer, 294, 174, 258, 1, 6, ncols, alpha);

            simdgeo::geom_p_y(buffer, 312, 174, 258, 1, 6, ncols, alpha);

            simdgeo::geom_p_z(buffer, 330, 174, 258, 1, 6, ncols, alpha);

            simdgeo::geom_p_x(buffer, 348, 162, 222, 1, 6, ncols, alpha);

            simdgeo::geom_p_y(buffer, 366, 162, 222, 1, 6, ncols, alpha);

            simdgeo::geom_p_z(buffer, 384, 162, 222, 1, 6, ncols, alpha);

            simdfunc::contract_primitives(buffer, 402, 348, 54, ncols);

            simdfunc::contract_primitives(buffer, 456, 294, 54, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 510, 456, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 510, 5, nmax);

    simdtrf::transform_d_inner(buffer, 510, 474, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 15 * nvalues, nvalues, buffer, 510, 5, nmax);

    simdtrf::transform_d_inner(buffer, 510, 492, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 30 * nvalues, nvalues, buffer, 510, 5, nmax);

    simdtrf::transform_d_inner(buffer, 510, 402, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 45 * nvalues, nvalues, buffer, 510, 5, nmax);

    simdtrf::transform_d_inner(buffer, 510, 420, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 60 * nvalues, nvalues, buffer, 510, 5, nmax);

    simdtrf::transform_d_inner(buffer, 510, 438, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 75 * nvalues, nvalues, buffer, 510, 5, nmax);
}

}  // namespace simdt2ceri
