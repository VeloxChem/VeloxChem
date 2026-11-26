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


#include "SimdElectronRepulsionGeom10RsRecSG.hpp"

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
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryS1.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_sg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_sg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 538, 448, 90, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5}, ncols,
                                                fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 12, {1, 2, 3, 4, 5}, ncols, fj,
                                            mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 24, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 27, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 30, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 33, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 36, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 39, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 42, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 45, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 48, 3, 10, 18, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 57, 3, 16, 21, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 66, 3, 7, 8, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 72, 3, 8, 9, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 78, 3, 9, 10, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 84, 3, 13, 14, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 90, 3, 14, 15, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 96, 3, 15, 16, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 102, 0, 3, 30, 78, 48, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 120, 0, 3, 42, 96, 57, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 138, 3, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 148, 3, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 158, 3, 36, 39, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 168, 3, 39, 42, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 178, 0, 3, 72, 148, 102, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 208, 0, 3, 90, 168, 120, ncols, p);

            compute_prim_sg_electron_repulsion_0(buffer, 238, 3, 66, 72, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 253, 3, 84, 90, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 268, 0, 3, 138, 238, 178, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 313, 0, 3, 158, 253, 208, ncols, p);

            simdgeo::geom_s_x(buffer, 358, 313, 1, 15, ncols, alpha);

            simdgeo::geom_s_y(buffer, 373, 313, 1, 15, ncols, alpha);

            simdgeo::geom_s_z(buffer, 388, 313, 1, 15, ncols, alpha);

            simdgeo::geom_s_x(buffer, 403, 268, 1, 15, ncols, alpha);

            simdgeo::geom_s_y(buffer, 418, 268, 1, 15, ncols, alpha);

            simdgeo::geom_s_z(buffer, 433, 268, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 448, 403, 45, ncols);

            simdfunc::contract_primitives(buffer, 493, 358, 45, ncols);
        }
    }

    simdtrf::transform_g_outer(values, nvalues, buffer, 493, 1, nmax);

    simdtrf::transform_g_outer(values + 9 * nvalues, nvalues, buffer, 508, 1, nmax);

    simdtrf::transform_g_outer(values + 18 * nvalues, nvalues, buffer, 523, 1, nmax);

    simdtrf::transform_g_outer(values + 27 * nvalues, nvalues, buffer, 448, 1, nmax);

    simdtrf::transform_g_outer(values + 36 * nvalues, nvalues, buffer, 463, 1, nmax);

    simdtrf::transform_g_outer(values + 45 * nvalues, nvalues, buffer, 478, 1, nmax);
}

}  // namespace simdt2ceri
