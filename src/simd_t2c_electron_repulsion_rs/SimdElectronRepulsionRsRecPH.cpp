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


#include "SimdElectronRepulsionRsRecPH.hpp"

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
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_ph_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ph_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 755, 596, 126, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 26, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 29, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 32, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 35, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 38, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 41, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 44, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 47, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 50, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 53, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 56, 3, 11, 20, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 65, 3, 18, 23, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 74, 3, 7, 8, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 80, 3, 8, 9, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 86, 3, 9, 10, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 92, 3, 10, 11, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 98, 3, 14, 15, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 104, 3, 15, 16, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 110, 3, 16, 17, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 116, 3, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 122, 0, 3, 35, 92, 56, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 140, 0, 3, 50, 116, 65, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 158, 3, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 168, 3, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 178, 3, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 188, 3, 41, 44, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 198, 3, 44, 47, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 208, 3, 47, 50, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 218, 0, 3, 86, 178, 122, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 248, 0, 3, 110, 208, 140, ncols, p);

            compute_prim_sg_electron_repulsion_0(buffer, 278, 3, 74, 80, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 293, 3, 80, 86, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 308, 3, 98, 104, 198, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 323, 3, 104, 110, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 338, 0, 3, 168, 293, 218, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 383, 0, 3, 198, 323, 248, ncols, p);

            compute_prim_sh_electron_repulsion_0(buffer, 428, 3, 158, 168, 293, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 449, 3, 188, 198, 323, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 470, 0, 3, 278, 428, 338, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 533, 0, 3, 308, 449, 383, ncols, p);

            simdfunc::contract_primitives(buffer, 596, 470, 126, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 722, 659, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 722, 11, nmax);

    simdtrf::transform_h_inner(buffer, 722, 596, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 33 * nvalues, nvalues, buffer, 722, 11, nmax);
}

}  // namespace simdt2ceri
