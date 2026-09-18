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


#include "SimdElectronRepulsionGeom10RsRecGP.hpp"

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
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_gp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_gp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1289, 974, 270, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 13, {1, 2, 3, 4, 5, 6}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 7, 8, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 8, 9, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 9, 10, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 10, 11, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 14, 15, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 15, 16, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 16, 17, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 17, 18, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 98, 0, 20, 23, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 108, 0, 23, 26, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 118, 0, 26, 29, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 128, 0, 35, 38, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 38, 41, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 41, 44, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 158, 0, 50, 56, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 173, 0, 56, 62, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 188, 0, 74, 80, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 203, 0, 80, 86, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 218, 0, 98, 108, 173, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 239, 0, 128, 138, 203, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 260, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 269, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 278, 3, 16, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 287, 3, 18, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 296, 0, 3, 23, 260, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 314, 0, 3, 29, 269, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 332, 0, 3, 38, 278, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 350, 0, 3, 44, 287, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 368, 0, 3, 50, 296, 98, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 398, 0, 3, 62, 314, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 428, 0, 3, 74, 332, 128, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 458, 0, 3, 86, 350, 148, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 488, 0, 3, 108, 398, 173, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 533, 0, 3, 138, 458, 203, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 578, 0, 3, 158, 488, 218, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 641, 0, 3, 188, 533, 239, ncols, p);

            simdgeo::geom_g_x(buffer, 704, 428, 641, 1, 3, ncols, alpha);

            simdgeo::geom_g_y(buffer, 749, 428, 641, 1, 3, ncols, alpha);

            simdgeo::geom_g_z(buffer, 794, 428, 641, 1, 3, ncols, alpha);

            simdgeo::geom_g_x(buffer, 839, 368, 578, 1, 3, ncols, alpha);

            simdgeo::geom_g_y(buffer, 884, 368, 578, 1, 3, ncols, alpha);

            simdgeo::geom_g_z(buffer, 929, 368, 578, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 974, 839, 135, ncols);

            simdfunc::contract_primitives(buffer, 1109, 704, 135, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 1244, 1109, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 1244, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1244, 1154, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 27 * nvalues, nvalues, buffer, 1244, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1244, 1199, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 54 * nvalues, nvalues, buffer, 1244, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1244, 974, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 81 * nvalues, nvalues, buffer, 1244, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1244, 1019, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 108 * nvalues, nvalues, buffer, 1244, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1244, 1064, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 135 * nvalues, nvalues, buffer, 1244, 3, nmax);
}

}  // namespace simdt2ceri
