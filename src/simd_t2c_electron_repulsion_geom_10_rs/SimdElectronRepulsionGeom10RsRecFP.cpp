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


#include "SimdElectronRepulsionGeom10RsRecFP.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryF1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_fp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_fp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 796, 586, 180, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 42, 0, 7, 8, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 8, 9, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 9, 10, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 13, 14, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 14, 15, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 15, 16, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 78, 0, 18, 21, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 88, 0, 21, 24, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 98, 0, 30, 33, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 108, 0, 33, 36, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 118, 0, 42, 48, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 133, 0, 60, 66, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_pp_electron_repulsion_0(buffer, 148, 3, 8, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 157, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 166, 3, 14, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 175, 3, 16, 39, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 184, 0, 3, 18, 148, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 202, 0, 3, 24, 157, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 220, 0, 3, 30, 166, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 238, 0, 3, 36, 175, 72, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 256, 0, 3, 48, 202, 88, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 286, 0, 3, 66, 238, 108, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 316, 0, 3, 78, 256, 118, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 361, 0, 3, 98, 286, 133, ncols, p);

            simdgeo::geom_f_x(buffer, 406, 220, 361, 1, 3, ncols, alpha);

            simdgeo::geom_f_y(buffer, 436, 220, 361, 1, 3, ncols, alpha);

            simdgeo::geom_f_z(buffer, 466, 220, 361, 1, 3, ncols, alpha);

            simdgeo::geom_f_x(buffer, 496, 184, 316, 1, 3, ncols, alpha);

            simdgeo::geom_f_y(buffer, 526, 184, 316, 1, 3, ncols, alpha);

            simdgeo::geom_f_z(buffer, 556, 184, 316, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 586, 496, 90, ncols);

            simdfunc::contract_primitives(buffer, 676, 406, 90, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 766, 676, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 766, 3, nmax);

    simdtrf::transform_p_inner(buffer, 766, 706, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 21 * nvalues, nvalues, buffer, 766, 3, nmax);

    simdtrf::transform_p_inner(buffer, 766, 736, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 42 * nvalues, nvalues, buffer, 766, 3, nmax);

    simdtrf::transform_p_inner(buffer, 766, 586, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 63 * nvalues, nvalues, buffer, 766, 3, nmax);

    simdtrf::transform_p_inner(buffer, 766, 616, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 84 * nvalues, nvalues, buffer, 766, 3, nmax);

    simdtrf::transform_p_inner(buffer, 766, 646, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 105 * nvalues, nvalues, buffer, 766, 3, nmax);
}

}  // namespace simdt2ceri
