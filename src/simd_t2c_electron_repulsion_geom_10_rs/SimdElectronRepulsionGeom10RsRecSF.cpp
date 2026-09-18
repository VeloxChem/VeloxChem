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


#include "SimdElectronRepulsionGeom10RsRecSF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryS1.hpp"
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_sf_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_sf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 318, 258, 60, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 22, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 25, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 28, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 31, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 34, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 37, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 40, 3, 9, 16, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 49, 3, 14, 19, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 58, 3, 7, 8, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 64, 3, 8, 9, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 70, 3, 12, 13, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 76, 3, 13, 14, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 82, 0, 3, 25, 64, 40, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 100, 0, 3, 34, 76, 49, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 118, 3, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 128, 3, 31, 34, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 138, 0, 3, 58, 118, 82, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 168, 0, 3, 70, 128, 100, ncols, p);

            simdgeo::geom_s_x(buffer, 198, 168, 1, 10, ncols, alpha);

            simdgeo::geom_s_y(buffer, 208, 168, 1, 10, ncols, alpha);

            simdgeo::geom_s_z(buffer, 218, 168, 1, 10, ncols, alpha);

            simdgeo::geom_s_x(buffer, 228, 138, 1, 10, ncols, alpha);

            simdgeo::geom_s_y(buffer, 238, 138, 1, 10, ncols, alpha);

            simdgeo::geom_s_z(buffer, 248, 138, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 258, 228, 30, ncols);

            simdfunc::contract_primitives(buffer, 288, 198, 30, ncols);
        }
    }

    simdtrf::transform_f_outer(values, nvalues, buffer, 288, 1, nmax);

    simdtrf::transform_f_outer(values + 7 * nvalues, nvalues, buffer, 298, 1, nmax);

    simdtrf::transform_f_outer(values + 14 * nvalues, nvalues, buffer, 308, 1, nmax);

    simdtrf::transform_f_outer(values + 21 * nvalues, nvalues, buffer, 258, 1, nmax);

    simdtrf::transform_f_outer(values + 28 * nvalues, nvalues, buffer, 268, 1, nmax);

    simdtrf::transform_f_outer(values + 35 * nvalues, nvalues, buffer, 278, 1, nmax);
}

}  // namespace simdt2ceri
