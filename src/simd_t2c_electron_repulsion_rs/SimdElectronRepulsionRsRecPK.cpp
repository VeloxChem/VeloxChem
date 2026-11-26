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


#include "SimdElectronRepulsionRsRecPK.hpp"

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
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_pk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_pk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1649, 1388, 216, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 15, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 30, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 33, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 36, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 39, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 42, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 45, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 48, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 51, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 54, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 57, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 60, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 63, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 66, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 69, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 72, 3, 13, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 81, 3, 22, 27, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 90, 3, 7, 8, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 96, 3, 8, 9, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 102, 3, 9, 10, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 108, 3, 10, 11, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 114, 3, 11, 12, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 120, 3, 12, 13, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 126, 3, 16, 17, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 132, 3, 17, 18, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 138, 3, 18, 19, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 144, 3, 19, 20, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 150, 3, 20, 21, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 156, 3, 21, 22, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 162, 0, 3, 45, 120, 72, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 180, 0, 3, 66, 156, 81, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 198, 3, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 208, 3, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 218, 3, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 228, 3, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 238, 3, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 248, 3, 51, 54, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 258, 3, 54, 57, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 268, 3, 57, 60, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 278, 3, 60, 63, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 288, 3, 63, 66, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 298, 0, 3, 114, 238, 162, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 328, 0, 3, 150, 288, 180, ncols, p);

            compute_prim_sg_electron_repulsion_0(buffer, 358, 3, 90, 96, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 373, 3, 96, 102, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 388, 3, 102, 108, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 403, 3, 108, 114, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 418, 3, 126, 132, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 433, 3, 132, 138, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 448, 3, 138, 144, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 463, 3, 144, 150, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 478, 0, 3, 228, 403, 298, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 523, 0, 3, 278, 463, 328, ncols, p);

            compute_prim_sh_electron_repulsion_0(buffer, 568, 3, 198, 208, 373, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 589, 3, 208, 218, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 610, 3, 218, 228, 403, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 631, 3, 248, 258, 433, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 652, 3, 258, 268, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 673, 3, 268, 278, 463, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 694, 0, 3, 388, 610, 478, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 757, 0, 3, 448, 673, 523, ncols, p);

            compute_prim_si_electron_repulsion_0(buffer, 820, 3, 358, 373, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 848, 3, 373, 388, 610, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 876, 3, 418, 433, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 904, 3, 433, 448, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 932, 0, 3, 589, 848, 694, ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 1016, 0, 3, 652, 904, 757, ncols, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1100, 3, 568, 589, 848, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1136, 3, 631, 652, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 1172, 0, 3, 820, 1100, 932, ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 1280, 0, 3, 876, 1136, 1016, ncols, p);

            simdfunc::contract_primitives(buffer, 1388, 1172, 216, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 1604, 1496, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 1604, 15, nmax);

    simdtrf::transform_k_inner(buffer, 1604, 1388, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 45 * nvalues, nvalues, buffer, 1604, 15, nmax);
}

}  // namespace simdt2ceri
