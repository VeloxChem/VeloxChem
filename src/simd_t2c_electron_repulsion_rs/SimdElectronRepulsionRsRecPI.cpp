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


#include "SimdElectronRepulsionRsRecPI.hpp"

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
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_pi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_pi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1139, 932, 168, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 28, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 31, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 34, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 37, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 40, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 43, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 46, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 49, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 52, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 55, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 58, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 61, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 64, 3, 12, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 73, 3, 20, 25, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 82, 3, 7, 8, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 88, 3, 8, 9, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 94, 3, 9, 10, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 100, 3, 10, 11, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 106, 3, 11, 12, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 112, 3, 15, 16, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 118, 3, 16, 17, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 124, 3, 17, 18, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 130, 3, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 136, 3, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 142, 0, 3, 40, 106, 64, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 160, 0, 3, 58, 136, 73, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 178, 3, 28, 31, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 188, 3, 31, 34, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 198, 3, 34, 37, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 208, 3, 37, 40, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 218, 3, 46, 49, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 228, 3, 49, 52, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 238, 3, 52, 55, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 248, 3, 55, 58, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 258, 0, 3, 100, 208, 142, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 288, 0, 3, 130, 248, 160, ncols, p);

            compute_prim_sg_electron_repulsion_0(buffer, 318, 3, 82, 88, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 333, 3, 88, 94, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 348, 3, 94, 100, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 363, 3, 112, 118, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 378, 3, 118, 124, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 393, 3, 124, 130, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 408, 0, 3, 198, 348, 258, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 453, 0, 3, 238, 393, 288, ncols, p);

            compute_prim_sh_electron_repulsion_0(buffer, 498, 3, 178, 188, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 519, 3, 188, 198, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 540, 3, 218, 228, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 561, 3, 228, 238, 393, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 582, 0, 3, 333, 519, 408, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 645, 0, 3, 378, 561, 453, ncols, p);

            compute_prim_si_electron_repulsion_0(buffer, 708, 3, 318, 333, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 736, 3, 363, 378, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 764, 0, 3, 498, 708, 582, ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 848, 0, 3, 540, 736, 645, ncols, p);

            simdfunc::contract_primitives(buffer, 932, 764, 168, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 1100, 1016, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 1100, 13, nmax);

    simdtrf::transform_i_inner(buffer, 1100, 932, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 39 * nvalues, nvalues, buffer, 1100, 13, nmax);
}

}  // namespace simdt2ceri
