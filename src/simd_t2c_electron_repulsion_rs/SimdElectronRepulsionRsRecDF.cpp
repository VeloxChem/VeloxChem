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


#include "SimdElectronRepulsionRsRecDF.hpp"

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
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_df_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_df_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 788, 626, 120, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 78, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 81, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 84, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 87, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 90, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 93, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 96, 3, 8, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 105, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 114, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 123, 3, 14, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 132, 3, 15, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 141, 3, 16, 39, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 150, 0, 3, 18, 96, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 168, 0, 3, 21, 105, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 186, 0, 3, 24, 114, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 204, 0, 3, 30, 123, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 222, 0, 3, 33, 132, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 240, 0, 3, 36, 141, 72, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 258, 3, 8, 9, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 264, 3, 9, 10, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 270, 3, 14, 15, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 276, 3, 15, 16, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 282, 0, 3, 78, 258, 105, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 300, 0, 3, 81, 264, 114, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 318, 0, 3, 87, 270, 132, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 336, 0, 3, 90, 276, 141, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 354, 0, 3, 105, 300, 42, 48, 186, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 390, 0, 3, 132, 336, 60, 66, 240, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 426, 3, 78, 81, 264, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 436, 3, 87, 90, 276, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 446, 0, 3, 258, 426, 300, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 476, 0, 3, 270, 436, 336, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 506, 0, 3, 282, 446, 150, 168, 354,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 566, 0, 3, 318, 476, 204, 222, 390,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 626, 506, 120, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 746, 686, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 746, 7, nmax);

    simdtrf::transform_f_inner(buffer, 746, 626, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 35 * nvalues, nvalues, buffer, 746, 7, nmax);
}

}  // namespace simdt2ceri
