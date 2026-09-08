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


#include "SimdElectronRepulsionRecDG.hpp"

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
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_dg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 775, 631, 90, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 29, 0, 7, 8, 14, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 35, 0, 8, 9, 17, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 41, 0, 9, 10, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 47, 0, 10, 11, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 11, 12, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 59, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 62, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 65, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 68, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 71, 3, 9, 17, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 80, 3, 10, 20, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 89, 3, 11, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 98, 3, 12, 26, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 107, 0, 3, 17, 80, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 125, 0, 3, 20, 89, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 143, 0, 3, 23, 98, 53, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 161, 3, 9, 10, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 167, 3, 10, 11, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 173, 3, 11, 12, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 179, 0, 3, 59, 161, 80, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 197, 0, 3, 62, 167, 89, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 215, 0, 3, 65, 173, 98, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 233, 0, 3, 71, 179, 29, 35, 107, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 269, 0, 3, 80, 197, 35, 41, 125, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 305, 0, 3, 89, 215, 41, 47, 143, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 341, 3, 59, 62, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 351, 3, 62, 65, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 361, 0, 3, 161, 341, 197, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 391, 0, 3, 167, 351, 215, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 421, 0, 3, 197, 391, 107, 125, 305,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 481, 3, 161, 167, 351, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 496, 0, 3, 341, 481, 391, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 541, 0, 3, 361, 496, 233, 269, 421,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 631, 541, 90, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 721, 631, 6, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 721, 9, nmax);
}

}  // namespace simdt2ceri
